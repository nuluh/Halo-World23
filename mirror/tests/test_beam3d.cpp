// Validation tests for the Beam3D element against closed-form solutions
// (biaxial Euler-Bernoulli bending, St. Venant torsion, and a right-angle
// grillage with coupled bending/torsion). Statics are solved by dynamic
// relaxation: damped explicit stepping until the transient dies, mirroring
// test_beam2d.cpp (same checkNear pattern and relax() recipe).
//
// Build: see scripts/build_native_tests.sh
#include <openvfife/solver/struct_system.hpp>
#include <openvfife/element/frame.hpp>
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

using openvfife::StdArray3d;
using openvfife::StdArray6d;
using openvfife::DOF;
using openvfife::Particle;
using openvfife::BaseMaterial;
using openvfife::LinearElastic;
using openvfife::CustomSectionParas;
using openvfife::CustomSection;
using openvfife::Beam3D;
using openvfife::StructSystem;
using namespace std;

static int g_checks = 0;
static int g_failures = 0;

static Particle* addParticle(StructSystem& system, int id, double x, double y, double z)
{
    auto particle = std::make_unique<Particle>(id, x, y, z);
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
// read back from the particles/elements afterwards (identical recipe to
// test_beam2d.cpp's relax())
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

// common material for every case: RHO=7850, E=200e9, and G set so that
// G = E/2.6 (i.e. nu=0.3, since BaseMaterial::setNu computes
// G_ = E_/2/(1+nu_)).
// FLAG: setNu() reads the material's CURRENT E_ to compute G_, so setE()
// MUST be called before setNu(), otherwise G would be derived from E_=0.
static void makeSteel(LinearElastic &mat)
{
    mat.setE(2.0e11);
    mat.setDensity(7850.0);
    mat.setNu(0.3);   // must follow setE(): G_ = E_ / 2 / (1 + 0.3) = E/2.6
}

// ---------------------------------------------------------------- case A
// biaxial cantilever along global x: two independent bending planes
// (Izz -> global y, Iyy -> global z) driven one at a time, 8x Beam3D
static void caseBiaxialCantilever()
{
    printf("case A: biaxial cantilever, Py and Pz tip loads (8x Beam3D)\n");

    const double L = 3.0, E = 2.0e11, A = 1.0e-2, Izz = 8.0e-5, Iyy = 2.0e-5,
                 J = 5.0e-5, P = 1.0e4;
    const int ne = 8;
    // Legacy bepcrd third-point direction +Y preserves ey=+Y, ez=+Z.
    const StdArray3d ref = {0.0, 1.0, 0.0};

    // --- sub-case 1: tip load along global y, checked against Izz -------
    {
        StructSystem sys(1);
        vector<Particle*> ps;
        for (int i = 0; i <= ne; i++)
        {
            ps.push_back(addParticle(sys, i + 1, L * i / ne, 0, 0));
        }

        LinearElastic mat(1);
        makeSteel(mat);
        CustomSectionParas sp = {A, 0, 0, 0, 0, 0, 0, Iyy, Izz, 0, J};
        CustomSection sect(1, sp);

        for (int i = 0; i < ne; i++)
        {
            sys.addElement(std::make_unique<Beam3D>(
                i + 1, *ps[i], *ps[i + 1], mat, sect, ref));
        }

        DOF fixAll = {{true, true, true, true, true, true},
                      {0, 0, 0, 0, 0, 0}};
        sys.addConstraint(1, fixAll);

        double h = sys.autoTimeStep();
        StdArray6d f = {0, P, 0, 0, 0, 0};
        relax(sys, {{ne + 1, f}}, h, 160.0, 1.0);

        const StdArray6d& ut = ps[ne]->displace();
        checkNear("Py: tip Uy", ut[1], P * L * L * L / (3 * E * Izz),
                  2e-3);
        checkNear("Py: tip Rotz", ut[5], P * L * L / (2 * E * Izz), 2e-3);

        const StdArray6d& fr = ps[0]->force();
        checkNear("Py: root reaction Fy", -fr[1], -P, 1e-3);
        checkNear("Py: root reaction |Mz|", fabs(-fr[5]), P * L, 1e-3);

    }

    // --- sub-case 2: fresh model, tip load along global z, checked
    // against Iyy ---------------------------------------------------
    {
        StructSystem sys(1);
        vector<Particle*> ps;
        for (int i = 0; i <= ne; i++)
        {
            ps.push_back(addParticle(sys, i + 1, L * i / ne, 0, 0));
        }

        LinearElastic mat(1);
        makeSteel(mat);
        CustomSectionParas sp = {A, 0, 0, 0, 0, 0, 0, Iyy, Izz, 0, J};
        CustomSection sect(1, sp);

        for (int i = 0; i < ne; i++)
        {
            sys.addElement(std::make_unique<Beam3D>(
                i + 1, *ps[i], *ps[i + 1], mat, sect, ref));
        }

        DOF fixAll = {{true, true, true, true, true, true},
                      {0, 0, 0, 0, 0, 0}};
        sys.addConstraint(1, fixAll);

        double h = sys.autoTimeStep();
        StdArray6d f = {0, 0, P, 0, 0, 0};
        relax(sys, {{ne + 1, f}}, h, 160.0, 1.0);

        const StdArray6d& ut = ps[ne]->displace();
        checkNear("Pz: tip Uz", ut[2], P * L * L * L / (3 * E * Iyy),
                  2e-3);
        // the spec flags the sign of Roty as something to observe rather
        // than assert; the magnitude is checked and the sign printed.
        checkNear("Pz: tip |Roty|", fabs(ut[4]),
                  P * L * L / (2 * E * Iyy), 2e-3);
        printf("  (observed sign: tip Roty = % .6e)\n", ut[4]);

        const StdArray6d& fr = ps[0]->force();
        checkNear("Pz: root reaction Fz", -fr[2], -P, 1e-3);
        checkNear("Pz: root reaction |My|", fabs(-fr[4]), P * L, 1e-3);

    }
}

// ---------------------------------------------------------------- case B
// torsion rod: tip moment about x on the same cantilever geometry
static void caseTorsionRod()
{
    printf("case B: torsion rod, tip moment about x (8x Beam3D)\n");

    const double L = 3.0, A = 1.0e-2, Izz = 8.0e-5, Iyy = 2.0e-5, J = 5.0e-5,
                 Mx = 1.0e3;
    const int ne = 8;
    // Legacy bepcrd third-point direction +Y preserves ey=+Y, ez=+Z.
    const StdArray3d ref = {0.0, 1.0, 0.0};

    StructSystem sys(1);
    vector<Particle*> ps;
    for (int i = 0; i <= ne; i++)
    {
        ps.push_back(addParticle(sys, i + 1, L * i / ne, 0, 0));
    }

    LinearElastic mat(1);
    makeSteel(mat);
    double G = mat.G();
    CustomSectionParas sp = {A, 0, 0, 0, 0, 0, 0, Iyy, Izz, 0, J};
    CustomSection sect(1, sp);

    for (int i = 0; i < ne; i++)
    {
        sys.addElement(std::make_unique<Beam3D>(
            i + 1, *ps[i], *ps[i + 1], mat, sect, ref));
    }

    DOF fixAll = {{true, true, true, true, true, true}, {0, 0, 0, 0, 0, 0}};
    sys.addConstraint(1, fixAll);

    double h = sys.autoTimeStep();
    StdArray6d f = {0, 0, 0, Mx, 0, 0};
    relax(sys, {{ne + 1, f}}, h, 160.0, 1.0);

    const StdArray6d& ut = ps[ne]->displace();
    checkNear("tip Rotx", ut[3], Mx * L / (G * J), 2e-3);

    const StdArray6d& fr = ps[0]->force();
    checkNear("root torque reaction", -fr[3], -Mx, 1e-3);

}

// ---------------------------------------------------------------- case C
// right-angle grillage: AB along x, BC along y, both fixed-jointed at B;
// a tip load Pz at the free corner C bends BC directly and both bends and
// twists AB, coupling bending and torsional stiffness in the same check
static void caseGrillage()
{
    printf("case C: right-angle grillage, Pz at free corner (8x8 Beam3D)\n");

    const double a = 2.0, b = 2.0, E = 2.0e11, A = 1.0e-2, I = 2.0e-5,
                 J = 4.0e-5, P = 1.0e3;
    const int ne = 8;
    // Preserve the original physical section axes under legacy bepcrd:
    // AB (+X): ey=+Y, ez=+Z; BC (+Y): ey=-X, ez=+Z.
    const StdArray3d refAB = {0.0, 1.0, 0.0};
    const StdArray3d refBC = {-1.0, 0.0, 0.0};

    StructSystem sys(1);

    // member AB: A(0,0,0) -> B(a,0,0)
    vector<Particle*> psAB;
    for (int i = 0; i <= ne; i++)
    {
        psAB.push_back(addParticle(sys, i + 1, a * i / ne, 0, 0));
    }

    // member BC: B(a,0,0) -> C(a,b,0); B is the shared joint (last node of
    // AB, first node of BC)
    vector<Particle*> psBC;
    psBC.push_back(psAB[ne]);
    for (int i = 1; i <= ne; i++)
    {
        psBC.push_back(addParticle(sys, ne + 1 + i, a, b * i / ne, 0));
    }

    LinearElastic mat(1);
    makeSteel(mat);
    double G = mat.G();
    CustomSectionParas sp = {A, 0, 0, 0, 0, 0, 0, I, I, 0, J};
    CustomSection sect(1, sp);

    int eid = 1;
    for (int i = 0; i < ne; i++)
    {
        sys.addElement(std::make_unique<Beam3D>(
            eid++, *psAB[i], *psAB[i + 1], mat, sect, refAB));
    }
    for (int i = 0; i < ne; i++)
    {
        sys.addElement(std::make_unique<Beam3D>(
            eid++, *psBC[i], *psBC[i + 1], mat, sect, refBC));
    }

    DOF fixAll = {{true, true, true, true, true, true}, {0, 0, 0, 0, 0, 0}};
    sys.addConstraint(1, fixAll);

    double h = sys.autoTimeStep();
    StdArray6d f = {0, 0, -P, 0, 0, 0};
    relax(sys, {{psBC[ne]->id(), f}}, h, 200.0, 2.0);

    const StdArray6d& uc = psBC[ne]->displace();
    double expected = -(P * (a * a * a + b * b * b) / (3 * E * I)
                        + P * a * b * b / (G * J));
    checkNear("tip Uz at C", uc[2], expected, 5e-3);

    const StdArray6d& fr = psAB[0]->force();
    checkNear("root reaction Fz at A", -fr[2], P, 1e-3);

}

int main()
{
    caseBiaxialCantilever();
    caseTorsionRod();
    caseGrillage();

    printf("\n%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
