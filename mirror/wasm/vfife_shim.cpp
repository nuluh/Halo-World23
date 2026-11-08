// C ABI shim wrapping StructSystem for the browser (Emscripten/WASM).
//
// The JavaScript worker builds a model through vf_add_* calls, runs the
// explicit solver in chunks with vf_run_steps, and reads flat result
// buffers refreshed by vf_collect. One model/session at a time; vf_reset
// tears everything down.
//
// Conventions exposed to JS:
//   nodes    vf_add_node places x, y in the x-y plane (z = 0);
//            vf_add_node3 places a full spatial node
//   elements vf_add_element: kind 0 = Beam2D, 1 = Link2DLD (planar);
//            vf_add_element3: kind 2 = Beam3D, 3 = Link3DLD (spatial),
//            with (rx, ry, rz) as the Beam3D reference vector — all
//            zero selects the default reference (0, 0, 1)
//   loads    vf_set_load(fx, fy, mz) / vf_set_load6(fx, fy, fz, mx, my,
//            mz) [N, N*m], one aggregate per node; gravity as
//            acceleration via vf_set_gravity(gx, gy) or
//            vf_set_gravity3(gx, gy, gz)
//   results  per node, stride 6: disp/vel slots [ux, uy, uz, rx, ry,
//            rz]; react = minus the accumulated force on each FIXED
//            dof, 0 on free dofs, same 6 slots
//            per element, stride 8: [N, V1, V2, T, Mz_i, Mz_j, My_i,
//            My_j].
//              2D elements fill [N, V, 0, 0, Mi, Mj, 0, 0]: N
//              (+tension), V (= dM/dx of sagging-positive BMD), Mi, Mj
//              (internal end moments, ccw+); the sagging BMD ordinates
//              at the ends are -Mi and +Mj.
//              Beam3D fills every slot from its end forces: V1/V2 shear
//              along local y/z, T torsion about local x, Mz_i/Mz_j the
//              plane-1 (Izz) end moments, My_i/My_j the plane-2 (Iyy)
//              end moments (about MINUS local y — see the Beam3D class
//              comment).
//              Links fill [N, 0, 0, 0, 0, 0, 0, 0].
//   fracture  capacity rupture checks on internal force magnitudes
//             (buckling is not modeled); checks and scheduled kills
//             only run when enabled via vf_set_fracture (the worker
//             enables it for dynamic runs only, because static
//             dynamic-relaxation transients overshoot ~2x and would
//             false-trigger).  Killed elements stop restoring forces
//             on the following step while their node masses remain
//             (debris keeps flying; there is no ground contact).
//             Beam2D checks max(|Mi|, |Mj|) against capM; Beam3D checks
//             max(|Mz_i|, |Mz_j|, |My_i|, |My_j|) (torsion is not
//             checked); |N| vs capN applies to every element kind.
#include <openvfife/solver/struct_system.hpp>
#include <openvfife/element/frame.hpp>
#include <openvfife/element/experimental.hpp>
#include <openvfife/material/uniaxial.hpp>
#include <openvfife/section/section.hpp>
#include <array>
#include <cmath>
#include <cstring>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <set>

using openvfife::StdArray3d;
using openvfife::StdArray6d;
using openvfife::DOF;
using openvfife::Mass;
using openvfife::Particle;
using openvfife::BaseMaterial;
using openvfife::LinearElastic;
using openvfife::UniBilinear;
using openvfife::BaseSection;
using openvfife::CustomSectionParas;
using openvfife::CustomSection;
using openvfife::ElementForce;
using openvfife::BaseElement;
using openvfife::StructElement;
using openvfife::Link2DLD;
using openvfife::Link3DLD;
using openvfife::Beam2D;
using openvfife::Beam3D;
using openvfife::StructSystem;

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>

#define VF_EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define VF_EXPORT
#endif

namespace
{

struct NodalLoad
{
    double fx, fy, fz;
    double mx, my, mz;
    double t0;   // window start; t1 <= t0 means always active
    double t1;   // window end (exclusive)
};

std::unique_ptr<StructSystem> g_sys;
std::map<int, Particle*> g_nodes;
std::map<int, BaseElement*> g_elems;   // insertion-ordered by id (map)
std::map<int, NodalLoad> g_loads;
std::map<int, DOF> g_supports;
std::array<double, 3> g_gravity{};
double g_dt = 0.0;
double g_zeta = 0.0;
double g_time = 0.0;
double g_load_scale = 1.0;
double g_ke_max = 0.0;      // probe phase: running kinetic-energy max
long g_peak_step = 0;       // probe phase: step index of that max
long g_auto_steps = 0;      // steps taken in auto (adaptive) mode
double g_auto_zeta = 0.0;   // 0 while probing, then pi / t_peak
bool g_started = false;
std::string g_error;

// element-kill and capacity-based fracture globals:
// capacities are rupture checks on internal force magnitudes (buckling
// is not modeled); checks and scheduled kills only run when fracture is
// enabled (the worker enables it for dynamic runs only, because static
// dynamic-relaxation transients overshoot ~2x and would false-trigger).
// Killed elements stop restoring forces on the following step while
// their node masses remain (debris keeps flying; there is no ground
// contact).
struct KillAt { int elem; double time; bool done; };
struct Capacity { int elem; double ncap; double mcap; };  // <=0 means no limit
std::vector<KillAt> g_kill_schedule;
std::vector<Capacity> g_capacities;
std::set<int> g_killed;                 // element ids already killed
std::vector<int> g_broken_ids;          // report buffers, kill order
std::vector<double> g_broken_times;
int g_fracture_on = 0;

// result buffers (flat, read by JS through heap views)
std::vector<int> g_buf_node_ids;
std::vector<double> g_buf_disp;      // 6 per node: ux, uy, uz, rx, ry, rz
std::vector<double> g_buf_vel;       // 6 per node
std::vector<double> g_buf_react;     // 6 per node (0 for free dofs)
std::vector<int> g_buf_elem_ids;
std::vector<double> g_buf_forces;    // 8 per element: N, V1, V2, T,
                                     // Mz_i, Mz_j, My_i, My_j (see the
                                     // header comment for the per-kind
                                     // fill)

void applyLoads()
{
    const bool grav = (g_gravity[0] != 0.0) || (g_gravity[1] != 0.0) ||
                      (g_gravity[2] != 0.0);
    for (auto& it : g_nodes)
    {
        Particle* p = it.second;
        double fx = 0.0, fy = 0.0, fz = 0.0;
        double mx = 0.0, my = 0.0, mz = 0.0;
        auto ld = g_loads.find(it.first);
        if (ld != g_loads.end())
        {
            const NodalLoad& nl = ld->second;
            if (nl.t1 <= nl.t0 || (g_time >= nl.t0 && g_time < nl.t1))
            {
                fx = nl.fx;
                fy = nl.fy;
                fz = nl.fz;
                mx = nl.mx;
                my = nl.my;
                mz = nl.mz;
            }
        }
        if (grav)
        {
            const Mass& m = p->mass();
            fx += m.lumped_mass[0] * g_gravity[0];
            fy += m.lumped_mass[1] * g_gravity[1];
            fz += m.lumped_mass[2] * g_gravity[2];
        }
        if (fx != 0.0 || fy != 0.0 || fz != 0.0 ||
            mx != 0.0 || my != 0.0 || mz != 0.0)
        {
            StdArray6d f = {g_load_scale * fx, g_load_scale * fy,
                            g_load_scale * fz, g_load_scale * mx,
                            g_load_scale * my, g_load_scale * mz};
            p->addForce(f);
        }
    }
}

// kill an element (idempotent: no-op if already killed or not found)
static void killNow(int id)
{
    if (g_killed.count(id)) return;
    auto it = g_elems.find(id);
    if (it == g_elems.end()) return;
    g_sys->killElement(id);
    g_killed.insert(id);
    g_broken_ids.push_back(id);
    g_broken_times.push_back(g_time);
}

// check scheduled-kill times and capacity-based rupture; called once
// per step after g_time has been advanced so the recorded break time
// equals the end-of-step time
static void checkFracture()
{
    if (!g_fracture_on) return;

    // scheduled kills
    for (auto& k : g_kill_schedule)
    {
        if (!k.done && g_time >= k.time)
        {
            killNow(k.elem);
            k.done = true;
        }
    }

    // capacity-based fracture
    for (auto& c : g_capacities)
    {
        if (g_killed.count(c.elem)) continue;
        auto it = g_elems.find(c.elem);
        if (it == g_elems.end()) continue;

        Beam2D* b2 = dynamic_cast<Beam2D*>(it->second);
        Beam3D* b3 = dynamic_cast<Beam3D*>(it->second);
        double N = 0.0, Mmax = 0.0;
        if (b2)
        {
            const ElementForce& fi = b2->end_force_i();
            const ElementForce& fj = b2->end_force_j();
            N = fi.fx;
            Mmax = std::max(fabs(fi.mz), fabs(fj.mz));
        }
        else if (b3)
        {
            const ElementForce& fi = b3->end_force_i();
            const ElementForce& fj = b3->end_force_j();
            N = fi.fx;
            Mmax = std::max(std::max(fabs(fi.mz), fabs(fj.mz)),
                            std::max(fabs(fi.my), fabs(fj.my)));
        }
        else
        {
            StructElement* s = dynamic_cast<StructElement*>(it->second);
            N = s ? s->element_force().fx : 0.0;
        }

        if ((c.ncap > 0.0 && fabs(N) > c.ncap) ||
            (c.mcap > 0.0 && Mmax > c.mcap))
        {
            killNow(c.elem);
        }
    }
}

bool fail(const std::string& msg)
{
    g_error = msg;
    return false;
}

double kineticEnergy()
{
    double ke = 0.0;
    for (auto& it : g_nodes)
    {
        Particle* p = it.second;
        const StdArray6d& v = p->velocity();
        const Mass& m = p->mass();
        for (int i = 0; i < 3; i++)
        {
            ke += 0.5 * m.lumped_mass[i] * v[i] * v[i];
        }

        Eigen::Vector3d omega(v[3], v[4], v[5]);
        if (m.rotation_mass_matrix.squaredNorm() > 0.0)
        {
            ke += 0.5 * omega.dot(m.rotation_mass_matrix * omega);
        }
        else
        {
            for (int i = 3; i < 6; i++)
            {
                ke += 0.5 * m.lumped_mass[i] * v[i] * v[i];
            }
        }
    }
    return ke;
}

void restartFromRest()
{
    StdArray6d zero = {0, 0, 0, 0, 0, 0};
    for (auto& it : g_nodes)
    {
        Particle* p = it.second;
        p->setPreviousDisplace(p->displace());
        p->setVelocity(zero);
        p->setPreviousVelocity(zero);
    }
}

// shared reference resolution/validation for vf_add_element and
// vf_add_element3; fills the out-pointers and returns true on success
bool resolveElementRefs(int id, int n1, int n2, int matId, int sectId,
                        Particle** p1, Particle** p2,
                        const BaseMaterial** mat, const BaseSection** sect)
{
    if (g_elems.count(id)) return fail("duplicate element id");
    auto i1 = g_nodes.find(n1);
    auto i2 = g_nodes.find(n2);
    if (i1 == g_nodes.end() || i2 == g_nodes.end())
    {
        return fail("element references unknown node");
    }
    const auto& mats = g_sys->materials();
    const auto& sects = g_sys->sections();
    auto mit = mats.find(matId);
    auto sit = sects.find(sectId);
    if (mit == mats.end()) return fail("unknown material");
    if (sit == sects.end()) return fail("unknown section");
    const auto& c1 = i1->second->coordinate();
    const auto& c2 = i2->second->coordinate();
    const double dx = c2[0] - c1[0];
    const double dy = c2[1] - c1[1];
    const double dz = c2[2] - c1[2];
    if (std::sqrt(dx * dx + dy * dy + dz * dz) < 1e-12)
    {
        return fail("zero-length element");
    }
    if (mit->second->density() <= 0.0)
    {
        return fail("material density must be positive (explicit solver "
                    "derives nodal mass from element mass)");
    }
    *p1 = i1->second;
    *p2 = i2->second;
    *mat = mit->second.get();
    *sect = sit->second.get();
    return true;
}

// record a support and register/replace its constraint on the system
void storeSupport(int nid, const DOF& d)
{
    if (g_supports.count(nid))
    {
        g_supports[nid] = d;
        g_sys->changeConstraint(nid, d);
    }
    else
    {
        g_supports[nid] = d;
        g_sys->addConstraint(nid, d);
    }
}

} // namespace

extern "C"
{

VF_EXPORT const char* vf_error()
{
    return g_error.c_str();
}

VF_EXPORT void vf_reset()
{
    g_sys.reset();
    g_nodes.clear();
    g_elems.clear();
    g_loads.clear();
    g_supports.clear();
    g_gravity.fill(0.0);
    g_dt = g_zeta = g_time = 0.0;
    g_ke_max = 0.0;
    g_peak_step = 0;
    g_auto_steps = 0;
    g_auto_zeta = 0.0;
    g_load_scale = 1.0;
    g_started = false;
    g_error.clear();
    g_kill_schedule.clear();
    g_capacities.clear();
    g_killed.clear();
    g_broken_ids.clear();
    g_broken_times.clear();
    g_fracture_on = 0;
}

VF_EXPORT int vf_create()
{
    vf_reset();
    g_sys = std::make_unique<StructSystem>(1);
    return 0;
}

// full spatial node; nothing is fixed here — constraints come from
// supports, and dof activation comes from the attached elements
VF_EXPORT int vf_add_node3(int id, double x, double y, double z)
{
    if (!g_sys) return -1;
    if (g_nodes.count(id)) { fail("duplicate node id"); return -1; }
    auto particle = std::make_unique<Particle>(id, x, y, z);
    Particle* p = particle.get();
    g_sys->addParticle(std::move(particle));
    g_nodes[id] = p;
    return 0;
}

// planar node in the x-y plane: the out-of-plane dofs (Uz, Rotx, Roty)
// stay inactive because no planar element kind ever activates them,
// exactly as before the 3D API existed
VF_EXPORT int vf_add_node(int id, double x, double y)
{
    return vf_add_node3(id, x, y, 0.0);
}

VF_EXPORT int vf_add_node_mass(int id, double m, double mrot)
{
    auto it = g_nodes.find(id);
    if (it == g_nodes.end()) { fail("node mass: unknown node"); return -1; }
    StdArray6d mm = {m, m, m, 0, 0, mrot};
    it->second->addLumpedMass(mm);
    return 0;
}

VF_EXPORT int vf_add_material(int id, double e, double density)
{
    if (!g_sys) return -1;
    auto mat = std::make_unique<LinearElastic>(id);
    mat->setE(e);
    mat->setDensity(density);
    g_sys->addMaterial(std::move(mat));
    return 0;
}

VF_EXPORT int vf_add_material_bilinear(int id, double e, double et,
                                       double fy, double m, double density)
{
    if (!g_sys) return -1;
    auto mat = std::make_unique<UniBilinear>(id, e, et, fy, m);
    mat->setDensity(density);
    g_sys->addMaterial(std::move(mat));
    return 0;
}

// set the shear modulus of an existing material (it defaults to 0, which
// zeroes Beam3D's torsional stiffness GJ/L). Beam3D reads G when the
// element is constructed, so call this before vf_add_element3.
VF_EXPORT int vf_set_material_g(int id, double g)
{
    if (!g_sys) return -1;
    const auto& mats = g_sys->materials();
    auto it = mats.find(id);
    if (it == mats.end()) { fail("unknown material"); return -1; }
    it->second->setG(g);
    return 0;
}

// full spatial section: j <= 0 selects the CustomSection polar fallback
// J = Iyy + Izz
VF_EXPORT int vf_add_section3(int id, double a, double iyy, double izz,
                              double j)
{
    if (!g_sys) return -1;
    CustomSectionParas sp = {a, 0, 0, 0, 0, 0, 0, iyy, izz, 0, j};
    auto sect = std::make_unique<CustomSection>(id, sp);
    g_sys->addSection(std::move(sect));
    return 0;
}

VF_EXPORT int vf_add_section(int id, double a, double izz)
{
    return vf_add_section3(id, a, izz, izz, 0.0);
}

// kind: 0 = Beam2D (frame member), 1 = Link2DLD (truss member,
// co-rotational large displacement); spatial kinds go through
// vf_add_element3
VF_EXPORT int vf_add_element(int id, int kind, int n1, int n2,
                             int matId, int sectId)
{
    if (!g_sys) return -1;
    if (kind != 0 && kind != 1)
    {
        fail("unknown element kind");
        return -1;
    }
    Particle* p1 = nullptr;
    Particle* p2 = nullptr;
    const BaseMaterial* mat = nullptr;
    const BaseSection* sect = nullptr;
    if (!resolveElementRefs(id, n1, n2, matId, sectId, &p1, &p2, &mat, &sect))
    {
        return -1;
    }
    std::unique_ptr<BaseElement> element;
    if (kind == 0)
    {
        element = std::make_unique<Beam2D>(id, *p1, *p2, *mat, *sect);
    }
    else
    {
        element = std::make_unique<Link2DLD>(id, *p1, *p2, *mat, *sect);
    }
    BaseElement* observer = element.get();
    g_sys->addElement(std::move(element));
    g_elems[id] = observer;
    return 0;
}

// kind: 2 = Beam3D (spatial frame member), 3 = Link3DLD (spatial truss
// member, co-rotational large displacement); planar kinds go through
// vf_add_element. (rx, ry, rz) is the Beam3D reference vector that seeds
// the section triad: all zero selects the default reference (0, 0, 1);
// it is ignored for kind 3.
VF_EXPORT int vf_add_element3(int id, int kind, int n1, int n2,
                              int matId, int sectId,
                              double rx, double ry, double rz)
{
    if (!g_sys) return -1;
    if (kind != 2 && kind != 3)
    {
        fail("vf_add_element3: element kind must be 2 (Beam3D) or "
             "3 (Link3DLD)");
        return -1;
    }
    Particle* p1 = nullptr;
    Particle* p2 = nullptr;
    const BaseMaterial* mat = nullptr;
    const BaseSection* sect = nullptr;
    if (!resolveElementRefs(id, n1, n2, matId, sectId, &p1, &p2, &mat, &sect))
    {
        return -1;
    }
    std::unique_ptr<BaseElement> element;
    if (kind == 2)
    {
        if (rx == 0.0 && ry == 0.0 && rz == 0.0)
        {
            element = std::make_unique<Beam3D>(id, *p1, *p2, *mat, *sect);
        }
        else
        {
            StdArray3d ref = {rx, ry, rz};
            element = std::make_unique<Beam3D>(id, *p1, *p2, *mat, *sect, ref);
        }
    }
    else
    {
        element = std::make_unique<Link3DLD>(id, *p1, *p2, *mat, *sect);
    }
    BaseElement* observer = element.get();
    g_sys->addElement(std::move(element));
    g_elems[id] = observer;
    return 0;
}

// planar support: the out-of-plane dofs (Uz, Rotx, Roty) are always
// fixed, exactly as before the 3D API existed
VF_EXPORT int vf_add_support(int nid, int fixUx, int fixUy, int fixRz,
                             double ux, double uy, double rz)
{
    if (!g_sys) return -1;
    if (!g_nodes.count(nid)) { fail("support on unknown node"); return -1; }
    DOF d = {{fixUx != 0, fixUy != 0, true, true, true, fixRz != 0},
             {ux, uy, 0, 0, 0, rz}};
    storeSupport(nid, d);
    return 0;
}

// full spatial support: each dof is fixed (to zero) where its flag is
// nonzero; unlike vf_add_support nothing is fixed implicitly
VF_EXPORT int vf_add_support6(int nid, int fUx, int fUy, int fUz,
                              int fRx, int fRy, int fRz)
{
    if (!g_sys) return -1;
    if (!g_nodes.count(nid)) { fail("support on unknown node"); return -1; }
    DOF d = {{fUx != 0, fUy != 0, fUz != 0, fRx != 0, fRy != 0, fRz != 0},
             {0, 0, 0, 0, 0, 0}};
    storeSupport(nid, d);
    return 0;
}

VF_EXPORT int vf_set_load6(int nid, double fx, double fy, double fz,
                           double mx, double my, double mz)
{
    if (!g_nodes.count(nid)) { fail("load on unknown node"); return -1; }
    g_loads[nid] = {fx, fy, fz, mx, my, mz, 0.0, 0.0};
    return 0;
}

VF_EXPORT int vf_set_load(int nid, double fx, double fy, double mz)
{
    return vf_set_load6(nid, fx, fy, 0.0, 0.0, 0.0, mz);
}

// time-windowed load: active while t0 <= t < t1 (an impulse); pass
// t1 <= t0 for an always-on load. Only meaningful in dynamic runs. The
// window lives in the same per-node aggregate as vf_set_load6, so the
// windowing mechanism applies per node no matter which setter wrote
// the load last.
VF_EXPORT int vf_set_load_window(int nid, double fx, double fy, double mz,
                                 double t0, double t1)
{
    if (!g_nodes.count(nid)) { fail("load on unknown node"); return -1; }
    g_loads[nid] = {fx, fy, 0.0, 0.0, 0.0, mz, t0, t1};
    return 0;
}

VF_EXPORT void vf_set_gravity3(double gx, double gy, double gz)
{
    g_gravity[0] = gx;
    g_gravity[1] = gy;
    g_gravity[2] = gz;
}

VF_EXPORT void vf_set_gravity(double gx, double gy)
{
    vf_set_gravity3(gx, gy, 0.0);
}

VF_EXPORT void vf_set_load_scale(double s)
{
    g_load_scale = s;
}

// capacities are rupture checks on internal force magnitudes (buckling
// is not modeled); checks and scheduled kills only run when fracture is
// enabled (the worker enables it for dynamic runs only, because static
// dynamic-relaxation transients overshoot ~2x and would false-trigger).
// Killed elements stop restoring forces on the following step while
// their node masses remain (debris keeps flying; no ground contact).
// The moment check uses max(|Mi|, |Mj|) for Beam2D and
// max(|Mz_i|, |Mz_j|, |My_i|, |My_j|) for Beam3D (torsion not checked).
VF_EXPORT int vf_set_capacity(int elemId, double ncap, double mcap)
{
    if (!g_elems.count(elemId))
    {
        fail("vf_set_capacity: unknown element id " + std::to_string(elemId));
        return -1;
    }
    for (auto& c : g_capacities)
    {
        if (c.elem == elemId)
        {
            c.ncap = ncap;
            c.mcap = mcap;
            return 0;
        }
    }
    g_capacities.push_back({elemId, ncap, mcap});
    return 0;
}

VF_EXPORT int vf_schedule_kill(int elemId, double time)
{
    if (!g_elems.count(elemId))
    {
        fail("vf_schedule_kill: unknown element id " + std::to_string(elemId));
        return -1;
    }
    for (auto& k : g_kill_schedule)
    {
        if (k.elem == elemId)
        {
            k.time = time;
            k.done = false;
            return 0;
        }
    }
    g_kill_schedule.push_back({elemId, time, false});
    return 0;
}

VF_EXPORT void vf_set_fracture(int enabled)
{
    g_fracture_on = enabled;
}

VF_EXPORT int vf_broken_count()
{
    return static_cast<int>(g_broken_ids.size());
}

VF_EXPORT int* vf_broken_ids_ptr()
{
    return g_broken_ids.data();
}

VF_EXPORT double* vf_broken_times_ptr()
{
    return g_broken_times.data();
}

VF_EXPORT double vf_critical_dt()
{
    // generic over every element kind: StructSystem::addElement folds
    // each element's virtual calcCriticalTimeStep() (Beam3D included)
    // into the running minimum
    if (!g_sys) return 0.0;
    return g_sys->critical_time_step();
}

// validates the model and performs the first (starting) step
VF_EXPORT int vf_setup(double dt, double zeta)
{
    if (!g_sys) { fail("no system"); return -1; }
    if (g_nodes.empty()) { fail("model has no nodes"); return -1; }
    if (g_elems.empty()) { fail("model has no elements"); return -1; }
    if (dt <= 0.0) { fail("dt must be positive"); return -1; }

    // every active dof must carry mass, otherwise the central difference
    // divides by zero
    for (auto& it : g_nodes)
    {
        const auto& dofs = it.second->dof_index();
        const auto& m = it.second->mass();
        for (int i = 0; i < 3; i++)
        {
            if (dofs[i] && m.lumped_mass[i] <= 0.0)
            {
                fail("active translational dof with zero mass at node " +
                     std::to_string(it.first) +
                     " (check material density / section area)");
                return -1;
            }
        }

        const bool any_rotation = dofs[3] || dofs[4] || dofs[5];
        if (any_rotation)
        {
            const Eigen::Matrix3d& inertia = m.rotation_mass_matrix;
            const bool matrix_inertia = inertia.squaredNorm() > 0.0;
            if (matrix_inertia)
            {
                // Particle::solve() inverts the complete coupled rotational
                // inertia tensor.  Require it to be symmetric positive
                // definite instead of demanding obsolete scalar lumped
                // masses on Rotx/Roty/Rotz.
                const double minor1 = inertia(0, 0);
                const double minor2 =
                    inertia(0, 0) * inertia(1, 1) -
                    inertia(0, 1) * inertia(1, 0);
                const double minor3 = inertia.determinant();
                const double symmetry_error =
                    (inertia - inertia.transpose()).cwiseAbs().maxCoeff();
                if (symmetry_error > 1e-12 ||
                    minor1 <= 0.0 || minor2 <= 0.0 || minor3 <= 0.0)
                {
                    fail("active rotational dof with invalid inertia matrix at node " +
                         std::to_string(it.first));
                    return -1;
                }
            }
            else
            {
                for (int i = 3; i < 6; i++)
                {
                    if (dofs[i] && m.lumped_mass[i] <= 0.0)
                    {
                        fail("active rotational dof with zero inertia at node " +
                             std::to_string(it.first));
                        return -1;
                    }
                }
            }
        }
        // a moment load on a node without an active rotational dof
        // would be silently ignored
        auto ld = g_loads.find(it.first);
        if (ld != g_loads.end())
        {
            const NodalLoad& nl = ld->second;
            if ((nl.mx != 0.0 && !dofs[3]) ||
                (nl.my != 0.0 && !dofs[4]) ||
                (nl.mz != 0.0 && !dofs[5]))
            {
                fail("moment load on node " + std::to_string(it.first) +
                     " which has no rotational dof (truss-only node)");
                return -1;
            }
        }
    }

    g_dt = dt;
    g_zeta = zeta;
    g_time = 0.0;
    g_ke_max = 0.0;
    g_peak_step = 0;
    g_auto_steps = 0;
    g_auto_zeta = 0.0;
    g_killed.clear();
    g_broken_ids.clear();
    g_broken_times.clear();
    for (auto& k : g_kill_schedule) k.done = false;

    applyLoads();
    g_sys->solve(dt, zeta, true);
    g_sys->setInternalForce();
    g_started = true;
    return 0;
}

// advance n steps; returns n, or -1 if the solution went non-finite
// (time step too large / unstable model)
VF_EXPORT int vf_run_steps(int n)
{
    if (!g_sys || !g_started) { fail("vf_setup not called"); return -1; }
    for (int i = 0; i < n; i++)
    {
        g_sys->solve(g_dt, g_zeta);
        g_sys->clearParticleForce();
        applyLoads();
        g_sys->setInternalForce();
        g_time += g_dt;
        checkFracture();
    }
    // stability probe over all six dofs per node (inactive dofs stay 0)
    for (auto& it : g_nodes)
    {
        const StdArray6d& u = it.second->displace();
        for (int i = 0; i < 6; i++)
        {
            if (!std::isfinite(u[i]))
            {
                fail("solution diverged (reduce dt or increase damping)");
                return -1;
            }
        }
    }
    return n;
}

// zero velocities and restart from rest: resets the dynamic state while
// preserving displacements and accumulated internal forces
VF_EXPORT void vf_zero_velocities()
{
    if (!g_sys) return;
    restartFromRest();
    g_ke_max = 0.0;
    g_peak_step = 0;
    g_auto_steps = 0;
}

// Adaptive dynamic relaxation: undamped probe run to measure the dominant
// period from the first kinetic-energy peak (≈ quarter period under a
// suddenly applied load), then switch permanently to near-critical
// mass-proportional viscous damping at zeta = pi / t_peak.  The damped
// stepping uses Particle::solve(h, zeta) which is unconditionally stable
// for mass-proportional damping, so the only stability constraint is
// plain central-difference (z < 2), and the caller can use the same dt as
// ordinary viscous stepping (0.7 * vf_critical_dt()).
//
// WHY THIS REPLACES KINETIC RESETS: a fixed-cadence multi-reset kinetic
// scheme is unstable because each velocity zeroing can amplify a mode by
// up to 1/cos(w*dt/2) when it lands at a bad phase.  Near equilibrium the
// kinetic-energy signal is dominated by high-frequency ripple, so resets
// fire at a fixed cadence and the CD-plus-velocity-zeroing map amplifies
// whichever mode lands at a bad phase — the solution explodes.
VF_EXPORT int vf_run_steps_kinetic(int n)
{
    if (!g_sys || !g_started) { fail("vf_setup not called"); return -1; }
    for (int i = 0; i < n; i++)
    {
        g_sys->solve(g_dt, g_auto_zeta);
        g_sys->clearParticleForce();
        applyLoads();
        g_sys->setInternalForce();
        g_time += g_dt;
        checkFracture();
        ++g_auto_steps;

        if (g_auto_zeta == 0.0)
        {
            const double ke = kineticEnergy();
            if (ke > g_ke_max)
            {
                g_ke_max = ke;
                g_peak_step = g_auto_steps;
            }
            else if (g_ke_max > 0.0 && ke < 0.8 * g_ke_max &&
                     g_peak_step >= 20)
            {
                const double t_peak = g_peak_step * g_dt;
                g_auto_zeta = M_PI / t_peak;
                restartFromRest();
            }
        }
    }
    for (auto& it : g_nodes)
    {
        const StdArray6d& u = it.second->displace();
        for (int i = 0; i < 6; i++)
        {
            if (!std::isfinite(u[i]))
            {
                fail("solution diverged (reduce dt)");
                return -1;
            }
        }
    }
    return n;
}

VF_EXPORT double vf_time()
{
    return g_time;
}

VF_EXPORT int vf_node_count()
{
    return static_cast<int>(g_nodes.size());
}

VF_EXPORT int vf_elem_count()
{
    return static_cast<int>(g_elems.size());
}

// refresh all result buffers from the current solver state
VF_EXPORT void vf_collect()
{
    const size_t nn = g_nodes.size();
    const size_t ne = g_elems.size();
    g_buf_node_ids.resize(nn);
    g_buf_disp.resize(6 * nn);
    g_buf_vel.resize(6 * nn);
    g_buf_react.resize(6 * nn);
    g_buf_elem_ids.resize(ne);
    g_buf_forces.resize(8 * ne);

    size_t k = 0;
    for (auto& it : g_nodes)
    {
        Particle* p = it.second;
        const StdArray6d& u = p->displace();
        const StdArray6d& v = p->velocity();
        const StdArray6d& f = p->force();
        g_buf_node_ids[k] = it.first;
        const bool sup = g_supports.count(it.first) > 0;
        const DOF* d = sup ? &g_supports[it.first] : nullptr;
        for (int i = 0; i < 6; i++)
        {
            g_buf_disp[6 * k + i] = u[i];
            g_buf_vel[6 * k + i] = v[i];
            // on a constrained dof the accumulated force is minus the
            // support reaction; free dofs report 0 here (their residual
            // is available through vf_max_residual_*)
            g_buf_react[6 * k + i] = (d && d->key[i]) ? -f[i] : 0.0;
        }
        ++k;
    }

    k = 0;
    for (auto& it : g_elems)
    {
        g_buf_elem_ids[k] = it.first;
        if (g_killed.count(it.first))
        {
            for (int i = 0; i < 8; i++)
            {
                g_buf_forces[8 * k + i] = 0.0;
            }
            ++k;
            continue;
        }
        Beam2D* b2 = dynamic_cast<Beam2D*>(it.second);
        Beam3D* b3 = dynamic_cast<Beam3D*>(it.second);
        if (b2)
        {
            const ElementForce& fi = b2->end_force_i();
            const ElementForce& fj = b2->end_force_j();
            g_buf_forces[8 * k + 0] = fi.fx;
            g_buf_forces[8 * k + 1] = fi.sfy;
            g_buf_forces[8 * k + 2] = 0.0;
            g_buf_forces[8 * k + 3] = 0.0;
            g_buf_forces[8 * k + 4] = fi.mz;
            g_buf_forces[8 * k + 5] = fj.mz;
            g_buf_forces[8 * k + 6] = 0.0;
            g_buf_forces[8 * k + 7] = 0.0;
        }
        else if (b3)
        {
            const ElementForce& fi = b3->end_force_i();
            const ElementForce& fj = b3->end_force_j();
            g_buf_forces[8 * k + 0] = fi.fx;    // N
            g_buf_forces[8 * k + 1] = fi.sfy;   // V1 (plane 1, along ey)
            g_buf_forces[8 * k + 2] = fi.sfz;   // V2 (plane 2, along ez)
            g_buf_forces[8 * k + 3] = fi.mx;    // T
            g_buf_forces[8 * k + 4] = fi.mz;    // Mz_i
            g_buf_forces[8 * k + 5] = fj.mz;    // Mz_j
            g_buf_forces[8 * k + 6] = fi.my;    // My_i
            g_buf_forces[8 * k + 7] = fj.my;    // My_j
        }
        else
        {
            StructElement* s = dynamic_cast<StructElement*>(it.second);
            g_buf_forces[8 * k + 0] = s ? s->element_force().fx : 0.0;
            for (int i = 1; i < 8; i++)
            {
                g_buf_forces[8 * k + i] = 0.0;
            }
        }
        ++k;
    }
}

VF_EXPORT int* vf_node_ids_ptr() { return g_buf_node_ids.data(); }
VF_EXPORT double* vf_disp_ptr() { return g_buf_disp.data(); }
VF_EXPORT double* vf_vel_ptr() { return g_buf_vel.data(); }
VF_EXPORT double* vf_react_ptr() { return g_buf_react.data(); }
VF_EXPORT int* vf_elem_ids_ptr() { return g_buf_elem_ids.data(); }
VF_EXPORT double* vf_forces_ptr() { return g_buf_forces.data(); }

VF_EXPORT double vf_kinetic_energy()
{
    return kineticEnergy();
}

// out-of-balance force over free translational dofs (Ux, Uy, Uz)
VF_EXPORT double vf_max_residual_force()
{
    double r = 0.0;
    for (auto& it : g_nodes)
    {
        const StdArray6d& f = it.second->force();
        const auto& dofs = it.second->dof_index();
        const DOF* sup =
            g_supports.count(it.first) ? &g_supports[it.first] : nullptr;
        for (int i = 0; i < 3; i++)
        {
            if (dofs[i] && !(sup && sup->key[i]))
                r = std::max(r, std::fabs(f[i]));
        }
    }
    return r;
}

// out-of-balance moment over free rotational dofs (Rotx, Roty, Rotz)
VF_EXPORT double vf_max_residual_moment()
{
    double r = 0.0;
    for (auto& it : g_nodes)
    {
        const StdArray6d& f = it.second->force();
        const auto& dofs = it.second->dof_index();
        const DOF* sup =
            g_supports.count(it.first) ? &g_supports[it.first] : nullptr;
        for (int i = 3; i < 6; i++)
        {
            if (dofs[i] && !(sup && sup->key[i]))
                r = std::max(r, std::fabs(f[i]));
        }
    }
    return r;
}

} // extern "C"
