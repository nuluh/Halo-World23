// Validation tests for the Beam2D element (and the existing Link2D path)
// against closed-form Euler-Bernoulli / truss solutions. Statics are solved
// by dynamic relaxation: damped explicit stepping until the transient dies.
//
// Build: see scripts/build_native_tests.sh
#include <openvfife/solver/struct_system.hpp>
#include <openvfife/element/frame.hpp>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

using openvfife::StdArray6d;
using openvfife::DOF;
using openvfife::Particle;
using openvfife::LinearElastic;
using openvfife::CustomSectionParas;
using openvfife::CustomSection;
using openvfife::StructElement;
using openvfife::Link2D;
using openvfife::Beam2D;
using openvfife::StructSystem;
using namespace std;

static int g_checks = 0;
static int g_failures = 0;

static Particle* addParticle(StructSystem& system, int id, double x, double y)
{
    auto particle = std::make_unique<Particle>(id, x, y);
    Particle* observer = particle.get();
    system.addParticle(std::move(particle));
    return observer;
}

static void checkNear(const char* label, double got, double want, double tol)
{
    ++g_checks;
    double err = fabs(got - want);
    double ref = fabs(want) > 1e-12 ? fabs(want) : 1.0;
    bool ok = err <= tol * ref;
    if (!ok)
    {
        ++g_failures;
    }
    printf("  %-38s got % .6e  want % .6e  %s\n",
           label, got, want, ok ? "ok" : "FAIL");
}

// run a constant-load dynamic relaxation and return nothing; results are
// read back from the particles/elements afterwards
static void relax(StructSystem& sys, const vector<pair<int, StdArray6d>>& loads,
                  double h, double zeta, double endTime)
{
    int nStep = (int)ceil(endTime / h);
    for (int i = 0; i <= nStep; i++)
    {
        if (i == 0)
        {
            for (auto& ld : loads) sys.addExternalForce(ld.first, ld.second);
            sys.solve(h, zeta, true);
            sys.setInternalForce();
        }
        else
        {
            sys.solve(h, zeta);
            sys.clearParticleForce();
            for (auto& ld : loads) sys.addExternalForce(ld.first, ld.second);
            sys.setInternalForce();
        }
    }
}

// ---------------------------------------------------------------- case A
// axial bar of two Beam2D elements: N = P, tip u = PL/EA
static void caseAxialBar()
{
    printf("case A: axial bar (Beam2D)\n");
    StructSystem sys(1);

    double L = 10.0, E = 1.0e6, A = 1.0, I = 1.0e-2, rho = 10.0, P = 100.0;

    Particle* p1 = addParticle(sys, 1, 0, 0);
    Particle* p2 = addParticle(sys, 2, L / 2, 0);
    Particle* p3 = addParticle(sys, 3, L, 0);

    LinearElastic mat(1);
    mat.setE(E);
    mat.setDensity(rho);
    CustomSectionParas sp = {A, 0, 0, 0, 0, 0, 0, I, I, 0};
    CustomSection sect(1, sp);

    sys.addElement(std::make_unique<Beam2D>(1, *p1, *p2, mat, sect));
    sys.addElement(std::make_unique<Beam2D>(2, *p2, *p3, mat, sect));

    DOF fixAll = {{true, true, true, true, true, true}, {0, 0, 0, 0, 0, 0}};
    sys.addConstraint(1, fixAll);

    double h = sys.autoTimeStep();
    StdArray6d f = {P, 0, 0, 0, 0, 0};
    relax(sys, {{3, f}}, h, 30.0, 20.0);

    const StdArray6d& u3 = p3->displace();
    checkNear("tip axial displacement", u3[0], P * L / (E * A), 5e-3);

    const auto& elems = sys.elements();
    Beam2D* e1 = dynamic_cast<Beam2D*>(elems.at(1).get());
    checkNear("element 1 axial force N", e1->end_force_i().fx, P, 5e-3);

    const StdArray6d& fr = p1->force();
    checkNear("root reaction Fx (=-force)", -fr[0], -P, 5e-3);

}

// ---------------------------------------------------------------- case B
// cantilever, tip transverse load: d = PL^3/3EI, th = PL^2/2EI
static void caseCantileverTipLoad()
{
    printf("case B: cantilever with tip transverse load (8x Beam2D)\n");
    StructSystem sys(1);

    double L = 2.0, E = 1.0e8, b = 0.1, hgt = 0.1, rho = 10.0, P = 1.0;
    double A = b * hgt;
    double I = b * hgt * hgt * hgt / 12.0;
    int ne = 8;

    vector<Particle*> ps;
    for (int i = 0; i <= ne; i++)
    {
        ps.push_back(addParticle(sys, i + 1, L * i / ne, 0));
    }

    LinearElastic mat(1);
    mat.setE(E);
    mat.setDensity(rho);
    CustomSectionParas sp = {A, 0, 0, 0, 0, 0, 0, I, I, 0};
    CustomSection sect(1, sp);

    for (int i = 0; i < ne; i++)
    {
        sys.addElement(std::make_unique<Beam2D>(i + 1, *ps[i], *ps[i + 1], mat, sect));
    }

    DOF fixAll = {{true, true, true, true, true, true}, {0, 0, 0, 0, 0, 0}};
    sys.addConstraint(1, fixAll);

    double h = sys.autoTimeStep();
    StdArray6d f = {0, P, 0, 0, 0, 0};
    relax(sys, {{ne + 1, f}}, h, 160.0, 1.0);

    const StdArray6d& ut = ps[ne]->displace();
    checkNear("tip deflection Uy", ut[1], P * L * L * L / (3 * E * I), 1e-2);
    checkNear("tip rotation Rz", ut[5], P * L * L / (2 * E * I), 1e-2);

    const StdArray6d& fr = ps[0]->force();
    checkNear("root reaction Fy", -fr[1], -P, 1e-2);
    checkNear("root reaction Mz", -fr[5], -P * L, 1e-2);

    const auto& elems = sys.elements();
    Beam2D* e1 = dynamic_cast<Beam2D*>(elems.at(1).get());
    // sagging-positive bending moment at the fixed end is +P*L = -Mi
    checkNear("root bending moment (-Mi)", -e1->end_force_i().mz, P * L, 1e-2);
    checkNear("root shear V", e1->end_force_i().sfy, -P, 1e-2);

}

// ---------------------------------------------------------------- case C
// cantilever, tip moment: th = ML/EI, d = ML^2/2EI, constant M along beam
static void caseCantileverTipMoment()
{
    printf("case C: cantilever with tip moment (8x Beam2D)\n");
    StructSystem sys(1);

    double L = 2.0, E = 1.0e8, b = 0.1, hgt = 0.1, rho = 10.0, M = 1.0;
    double A = b * hgt;
    double I = b * hgt * hgt * hgt / 12.0;
    int ne = 8;

    vector<Particle*> ps;
    for (int i = 0; i <= ne; i++)
    {
        ps.push_back(addParticle(sys, i + 1, L * i / ne, 0));
    }

    LinearElastic mat(1);
    mat.setE(E);
    mat.setDensity(rho);
    CustomSectionParas sp = {A, 0, 0, 0, 0, 0, 0, I, I, 0};
    CustomSection sect(1, sp);

    for (int i = 0; i < ne; i++)
    {
        sys.addElement(std::make_unique<Beam2D>(i + 1, *ps[i], *ps[i + 1], mat, sect));
    }

    DOF fixAll = {{true, true, true, true, true, true}, {0, 0, 0, 0, 0, 0}};
    sys.addConstraint(1, fixAll);

    double h = sys.autoTimeStep();
    StdArray6d f = {0, 0, 0, 0, 0, M};
    relax(sys, {{ne + 1, f}}, h, 160.0, 1.0);

    const StdArray6d& ut = ps[ne]->displace();
    checkNear("tip rotation Rz", ut[5], M * L / (E * I), 1e-2);
    checkNear("tip deflection Uy", ut[1], M * L * L / (2 * E * I), 1e-2);

    const auto& elems = sys.elements();
    Beam2D* e4 = dynamic_cast<Beam2D*>(elems.at(4).get());
    // constant sagging moment +M along the member: -Mi = +Mj = M, V = 0
    checkNear("mid-span moment (-Mi)", -e4->end_force_i().mz, M, 1e-2);
    checkNear("mid-span moment (+Mj)", e4->end_force_j().mz, M, 1e-2);
    checkNear("mid-span shear V", e4->end_force_i().sfy, 0.0, 1e-2);

    const StdArray6d& fr = ps[0]->force();
    checkNear("root reaction Mz", -fr[5], -M, 1e-2);

}

// ---------------------------------------------------------------- case D
// simply supported beam, midspan point load: d = PL^3/48EI, M_mid = PL/4
static void caseSimplySupported()
{
    printf("case D: simply supported beam, midspan load (8x Beam2D)\n");
    StructSystem sys(1);

    double L = 4.0, E = 1.0e8, b = 0.1, hgt = 0.1, rho = 10.0, P = 1.0;
    double A = b * hgt;
    double I = b * hgt * hgt * hgt / 12.0;
    int ne = 8;

    vector<Particle*> ps;
    for (int i = 0; i <= ne; i++)
    {
        ps.push_back(addParticle(sys, i + 1, L * i / ne, 0));
    }

    LinearElastic mat(1);
    mat.setE(E);
    mat.setDensity(rho);
    CustomSectionParas sp = {A, 0, 0, 0, 0, 0, 0, I, I, 0};
    CustomSection sect(1, sp);

    for (int i = 0; i < ne; i++)
    {
        sys.addElement(std::make_unique<Beam2D>(i + 1, *ps[i], *ps[i + 1], mat, sect));
    }

    // left pin, right roller (Uz/Rotx/Roty are inactive in 2D anyway)
    DOF pin = {{true, true, true, true, true, false}, {0, 0, 0, 0, 0, 0}};
    DOF roller = {{false, true, true, true, true, false}, {0, 0, 0, 0, 0, 0}};
    sys.addConstraint(1, pin);
    sys.addConstraint(ne + 1, roller);

    double h = sys.autoTimeStep();
    StdArray6d f = {0, -P, 0, 0, 0, 0};
    relax(sys, {{ne / 2 + 1, f}}, h, 60.0, 2.0);

    const StdArray6d& um = ps[ne / 2]->displace();
    checkNear("midspan deflection Uy", um[1], -P * L * L * L / (48 * E * I),
              1e-2);

    const StdArray6d& fl = ps[0]->force();
    const StdArray6d& fright = ps[ne]->force();
    checkNear("left reaction Fy", -fl[1], P / 2, 1e-2);
    checkNear("right reaction Fy", -fright[1], P / 2, 1e-2);

    const auto& elems = sys.elements();
    Beam2D* e4 = dynamic_cast<Beam2D*>(elems.at(ne / 2).get());
    // downward midspan load sags the beam: M_mid = +PL/4 (= +Mj of elem 4)
    checkNear("midspan moment (+Mj)", e4->end_force_j().mz, P * L / 4, 1e-2);
    checkNear("left half shear V", e4->end_force_i().sfy, P / 2, 1e-2);

}

// ---------------------------------------------------------------- case E
// documented example 3: two-bar truss with Link2D, verifies upstream path
static void caseTrussExample3()
{
    printf("case E: two-bar truss (upstream Link2D, examples doc case 3)\n");
    StructSystem sys(1);

    Particle* p1 = addParticle(sys, 1, 0, 0);
    Particle* p2 = addParticle(sys, 2, 10, 0);
    Particle* p3 = addParticle(sys, 3, 10, 5);
    StdArray6d mass = {10.0, 10.0, 0, 0, 0, 0};
    p2->setLumpedMass(mass);

    LinearElastic mat1(1);
    mat1.setE(2.0e4);
    LinearElastic mat2(2);
    mat2.setE(1.0e4);

    CustomSectionParas sp = {1, 0, 0, 0, 0, 0, 0, 10, 10, 10};
    CustomSection sect(1, sp);

    sys.addElement(std::make_unique<Link2D>(1, *p1, *p2, mat1, sect));
    sys.addElement(std::make_unique<Link2D>(2, *p2, *p3, mat2, sect));

    DOF fixAll = {{true, true, true, true, true, true}, {0, 0, 0, 0, 0, 0}};
    sys.addConstraint(1, fixAll);
    sys.addConstraint(3, fixAll);

    // loads exactly as printed in examples/openVFIFE_examples.md example 3.
    // NOTE: the results table in that document does not correspond to these
    // load signs (it matches the opposite sign, and its reported reactions
    // are not in equilibrium with its own geometry); the expectations below
    // are hand statics for the coded loads:
    //   equilibrium of free node 2: N_a = N_b = -20 (compression),
    //   u2 = (N_a/k_a, -N_b/k_b) = (-0.01, +0.01) with k = EA/L = 2000,
    //   support 1: R = -(ext + bar) = (+40, 0),
    //   support 3: R = (-20, -20)
    StdArray6d f1 = {-20, 0, 0, 0, 0, 0};
    StdArray6d f2 = {-20, 20, 0, 0, 0, 0};
    StdArray6d f3 = {20, 0, 0, 0, 0, 0};
    relax(sys, {{1, f1}, {2, f2}, {3, f3}}, 0.01, 0.5, 50.0);

    const StdArray6d& u2 = p2->displace();
    checkNear("node 2 Ux", u2[0], -0.01, 1e-2);
    checkNear("node 2 Uy", u2[1], 0.01, 1e-2);

    const StdArray6d& r1 = p1->force();
    const StdArray6d& r3 = p3->force();
    checkNear("reaction 1 Fx", -r1[0], 40.0, 1e-2);
    checkNear("reaction 1 Fy", -r1[1], 0.0, 1e-2);
    checkNear("reaction 3 Fx", -r3[0], -20.0, 1e-2);
    checkNear("reaction 3 Fy", -r3[1], -20.0, 1e-2);

    const auto& elems = sys.elements();
    StructElement* e1 = dynamic_cast<StructElement*>(elems.at(1).get());
    StructElement* e2 = dynamic_cast<StructElement*>(elems.at(2).get());
    checkNear("bar a axial force", e1->element_force().fx, -20.0, 1e-2);
    checkNear("bar b axial force", e2->element_force().fx, -20.0, 1e-2);

}

int main()
{
    caseAxialBar();
    caseCantileverTipLoad();
    caseCantileverTipMoment();
    caseSimplySupported();
    caseTrussExample3();

    printf("\n%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
