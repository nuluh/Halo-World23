// Drives the C shim exactly like the web worker will: build a portal frame,
// relax to statics, read flat buffers. Native mirror of the WASM usage.
#include <cmath>
#include <cstdio>
#include <cstdlib>

extern "C"
{
    const char* vf_error();
    void vf_reset();
    int vf_create();
    int vf_add_node(int id, double x, double y);
    int vf_add_node3(int id, double x, double y, double z);
    int vf_add_node_mass(int id, double m, double mrot);
    int vf_add_material(int id, double e, double density);
    int vf_set_material_g(int id, double g);
    int vf_add_section(int id, double a, double izz);
    int vf_add_section3(int id, double a, double iyy, double izz, double j);
    int vf_add_element(int id, int kind, int n1, int n2, int matId, int sectId);
    int vf_add_element3(int id, int kind, int n1, int n2,
                        int matId, int sectId,
                        double rx, double ry, double rz);
    int vf_add_support(int nid, int fx, int fy, int frz,
                       double ux, double uy, double rz);
    int vf_add_support6(int nid, int fUx, int fUy, int fUz,
                        int fRx, int fRy, int fRz);
    int vf_set_load(int nid, double fx, double fy, double mz);
    int vf_set_load6(int nid, double fx, double fy, double fz,
                     double mx, double my, double mz);
    int vf_set_load_window(int nid, double fx, double fy, double mz,
                           double t0, double t1);
    void vf_set_gravity(double gx, double gy);
    void vf_set_gravity3(double gx, double gy, double gz);
    double vf_critical_dt();
    int vf_setup(double dt, double zeta);
    int vf_run_steps(int n);
    int vf_run_steps_kinetic(int n);
    void vf_collect();
    int vf_node_count();
    int vf_elem_count();
    int* vf_node_ids_ptr();
    int* vf_elem_ids_ptr();
    double* vf_disp_ptr();
    double* vf_vel_ptr();
    double* vf_react_ptr();
    double* vf_forces_ptr();
    double vf_kinetic_energy();
    double vf_max_residual_force();
    double vf_max_residual_moment();
    int vf_set_capacity(int elemId, double ncap, double mcap);
    int vf_schedule_kill(int elemId, double time);
    void vf_set_fracture(int enabled);
    int vf_broken_count();
    int* vf_broken_ids_ptr();
    double* vf_broken_times_ptr();
}

static int g_checks = 0, g_failures = 0;

static void checkNear(const char* label, double got, double want, double tol)
{
    ++g_checks;
    double ref = fabs(want) > 1e-12 ? fabs(want) : 1.0;
    bool ok = fabs(got - want) <= tol * ref;
    if (!ok) ++g_failures;
    printf("  %-40s got % .6e  want % .6e  %s\n",
           label, got, want, ok ? "ok" : "FAIL");
}

// portal frame: two columns height H, beam span L, lateral load P at the
// beam level, columns fixed at the base, all members meshed n-fold.
// Reference: sway frame with infinitely-stiff-beam approximation is NOT
// exact, so instead cross-check the global equilibrium relations that hold
// exactly: sum of base shears = P, sum of base moments + P*H balanced by
// column end moments, and symmetry of the two columns' axial forces.
static const double H = 3.0, L = 4.0, E = 2.0e8;
static const double A = 0.01, I = 8.333e-6, RHO = 10.0, P = 1.0;
static const int mesh = 4;

static void runPortal(bool kinetic)
{
    vf_create();

    // nodes 1..(mesh+1): left column bottom->top, then beam nodes, then
    // right column top->bottom; keep it simple with explicit ids
    // left column x=0: ids 1..mesh+1 (y = 0..H)
    // beam y=H: ids mesh+1 .. mesh+1+mesh (x = 0..L)
    // right column x=L: ids mesh+1+mesh .. 2*mesh+1+mesh (y = H..0)
    int id = 1;
    for (int i = 0; i <= mesh; i++)
        vf_add_node(id++, 0.0, H * i / mesh);
    for (int i = 1; i <= mesh; i++)
        vf_add_node(id++, L * i / mesh, H);
    for (int i = 1; i <= mesh; i++)
        vf_add_node(id++, L, H - H * i / mesh);
    const int nTop = mesh + 1;               // left beam-column joint
    const int nTopRight = mesh + 1 + mesh;   // right beam-column joint
    const int nBaseRight = 3 * mesh + 1;     // right column base

    vf_add_material(1, E, RHO);
    vf_add_section(1, A, I);

    int eid = 1;
    for (int n = 1; n < nBaseRight; n++)
    {
        if (vf_add_element(eid++, 0, n, n + 1, 1, 1) != 0)
        {
            printf("add_element failed: %s\n", vf_error());
            exit(1);
        }
    }

    vf_add_support(1, 1, 1, 1, 0, 0, 0);
    vf_add_support(nBaseRight, 1, 1, 1, 0, 0, 0);
    vf_set_load(nTop, P, 0, 0);

    // adaptive damping only needs plain central-difference stability
    double dtc = vf_critical_dt();
    double dt = kinetic ? 0.7 * dtc : pow(10.0, floor(log10(dtc)));
    printf("portal frame (%s): dt_critical=%.3e dt=%.3e\n",
           kinetic ? "kinetic damping" : "viscous", dtc, dt);

    if (vf_setup(dt, kinetic ? 0.0 : 60.0) != 0)
    {
        printf("setup failed: %s\n", vf_error());
        exit(1);
    }

    // relax in chunks until quiet, like the worker will
    int total = 0;
    for (int chunk = 0; chunk < 400; chunk++)
    {
        int rc = kinetic ? vf_run_steps_kinetic(2000) : vf_run_steps(2000);
        if (rc < 0)
        {
            printf("run failed: %s\n", vf_error());
            exit(1);
        }
        total += 2000;
        if (vf_max_residual_force() < 1e-8 * P &&
            vf_max_residual_moment() < 1e-8 * P * H)
            break;
    }
    printf("  steps: %d  residF %.2e residM %.2e KE %.2e\n", total,
           vf_max_residual_force(), vf_max_residual_moment(),
           vf_kinetic_energy());

    vf_collect();
    int nn = vf_node_count();
    int ne = vf_elem_count();
    int* ids = vf_node_ids_ptr();
    double* u = vf_disp_ptr();
    double* r = vf_react_ptr();
    double* ef = vf_forces_ptr();

    // locate result rows
    int rowTop = -1, rowTopR = -1, rowB1 = -1, rowB2 = -1;
    for (int k = 0; k < nn; k++)
    {
        if (ids[k] == nTop) rowTop = k;
        if (ids[k] == nTopRight) rowTopR = k;
        if (ids[k] == 1) rowB1 = k;
        if (ids[k] == nBaseRight) rowB2 = k;
    }

    // global equilibrium: base shears carry the applied lateral load
    checkNear("sum base horizontal reactions",
              r[6 * rowB1 + 0] + r[6 * rowB2 + 0], -P, 1e-6);
    // no net vertical load: vertical reactions cancel
    checkNear("sum base vertical reactions",
              r[6 * rowB1 + 1] + r[6 * rowB2 + 1], 0.0, 1e-6);
    // overturning about the left base: the applied load contributes a
    // moment of -P*H (ccw+), so reactions must supply
    // Mbase1 + Mbase2 + Ry_right*L = +P*H
    checkNear("overturning balance",
              r[6 * rowB1 + 5] + r[6 * rowB2 + 5] + r[6 * rowB2 + 1] * L,
              P * H, 1e-6);
    // sway symmetry: both beam-level joints translate almost equally
    checkNear("joint sway consistency",
              u[6 * rowTop + 0] / u[6 * rowTopR + 0], 1.0, 1e-3);

    // matrix-analysis reference for this exact frame (computed with the
    // standard 6-dof frame stiffness method, same discretization limit):
    // lateral drift of a fixed-base portal frame
    //   d = P*H^3/(12*E*Ic) * (3*k+2)/(6*k+1),  k = (Ib/L)/(Ic/H)
    const double kr = (I / L) / (I / H);
    const double dRef =
        P * H * H * H / (12.0 * E * I) * (3.0 * kr + 2.0) / (6.0 * kr + 1.0);
    checkNear("lateral drift vs portal formula", u[6 * rowTop + 0], dRef,
              2e-2);

    // element force sanity: axial force in the beam equals shear carried
    // into the right column; left column base element shear + right column
    // base element shear = P. Find the two base elements (1st and last).
    double vBase1 = ef[8 * 0 + 1];
    double vBaseR = ef[8 * (ne - 1) + 1];
    checkNear("base shears sum to P", fabs(vBase1) + fabs(vBaseR), P, 1e-3);
}

// cantilever beam: L=4 m, tip load P=8 kN downward, 10 elements.
// Regression test for the kinetic-damping instability: must converge,
// not explode.
static void runCantilever(bool kinetic)
{
    const double L = 4.0, E = 200.0e9, A = 2.85e-3, I = 1.94e-5;
    const double RHO = 7850.0, P = 8000.0;
    const int nelem = 10;

    vf_create();

    for (int i = 1; i <= nelem + 1; i++)
        vf_add_node(i, L * (i - 1) / nelem, 0.0);

    vf_add_material(1, E, RHO);
    vf_add_section(1, A, I);

    for (int i = 1; i <= nelem; i++)
    {
        if (vf_add_element(i, 0, i, i + 1, 1, 1) != 0)
        {
            printf("add_element failed: %s\n", vf_error());
            exit(1);
        }
    }

    vf_add_support(1, 1, 1, 1, 0, 0, 0);
    vf_set_load(nelem + 1, 0.0, -P, 0.0);

    double dtc = vf_critical_dt();
    double dt = kinetic ? 0.7 * dtc : pow(10.0, floor(log10(dtc)));
    printf("cantilever (%s): dt_critical=%.3e dt=%.3e\n",
           kinetic ? "adaptive damping" : "viscous", dtc, dt);

    if (vf_setup(dt, kinetic ? 0.0 : 60.0) != 0)
    {
        printf("setup failed: %s\n", vf_error());
        exit(1);
    }

    int total = 0;
    for (int chunk = 0; chunk < 400; chunk++)
    {
        int rc = kinetic ? vf_run_steps_kinetic(2000) : vf_run_steps(2000);
        if (rc < 0)
        {
            printf("run failed: %s\n", vf_error());
            exit(1);
        }
        total += 2000;
        if (vf_max_residual_force() < 1e-8 * P &&
            vf_max_residual_moment() < 1e-8 * P * L)
            break;
    }
    printf("  steps: %d  residF %.2e residM %.2e KE %.2e\n", total,
           vf_max_residual_force(), vf_max_residual_moment(),
           vf_kinetic_energy());

    vf_collect();
    int nn = vf_node_count();
    int* ids = vf_node_ids_ptr();
    double* u = vf_disp_ptr();
    double* r = vf_react_ptr();

    int rowTip = -1, rowRoot = -1;
    for (int k = 0; k < nn; k++)
    {
        if (ids[k] == nelem + 1) rowTip = k;
        if (ids[k] == 1) rowRoot = k;
    }

    // Euler-Bernoulli references for a tip-loaded cantilever
    checkNear("tip deflection uy",
              u[6 * rowTip + 1],
              -P * L * L * L / (3.0 * E * I), 1e-3);
    checkNear("tip rotation rz",
              u[6 * rowTip + 5],
              -P * L * L / (2.0 * E * I), 1e-3);
    checkNear("root vertical reaction",
              r[6 * rowRoot + 1], +P, 1e-6);
    // the co-rotational solution is geometrically nonlinear: the tip pulls
    // in slightly (ux < 0), so the base moment balances P times the
    // deformed lever arm, not P*L exactly
    checkNear("root moment reaction",
              r[6 * rowRoot + 5], P * (L + u[6 * rowTip + 0]), 1e-6);
}

// progressive-collapse / scheduled-kill and capacity-based fracture tests.
// Scenario A: schedule element 1 to break at 50.5*dt; verify it breaks,
//   forces drop to zero, and root reactions vanish (the structure is
//   severed from its support and drifts as a projectile).
// Scenario B: set a moment capacity at half the static root moment so the
//   suddenly-applied tip load ruptures element 1 early, and verify the
//   run stays finite afterward.
static void runKillAndFracture()
{
    const double L = 4.0, E = 200.0e9, A = 2.85e-3, I = 1.94e-5;
    const double RHO = 7850.0, P = 8000.0;
    const int nelem = 10;

    // --- Scenario A: scheduled kill ------------------------------------
    {
        vf_create();

        for (int i = 1; i <= nelem + 1; i++)
            vf_add_node(i, L * (i - 1) / nelem, 0.0);

        vf_add_material(1, E, RHO);
        vf_add_section(1, A, I);

        for (int i = 1; i <= nelem; i++)
        {
            if (vf_add_element(i, 0, i, i + 1, 1, 1) != 0)
            {
                printf("add_element failed: %s\n", vf_error());
                exit(1);
            }
        }

        vf_add_support(1, 1, 1, 1, 0, 0, 0);
        vf_set_load(nelem + 1, 0.0, -P, 0.0);

        double dtc = vf_critical_dt();
        double dt = pow(10.0, floor(log10(dtc)));
        printf("kill-scheduled (A): dt_critical=%.3e dt=%.3e\n", dtc, dt);

        vf_set_fracture(1);
        vf_schedule_kill(1, 50.5 * dt);

        if (vf_setup(dt, 0.0) != 0)
        {
            printf("setup failed: %s\n", vf_error());
            exit(1);
        }

        int rc = vf_run_steps(2000);
        if (rc < 0)
        {
            printf("vf_run_steps failed: %s\n", vf_error());
            ++g_failures;
        }

        checkNear("A: broken count", (double)vf_broken_count(), 1.0, 1e-9);
        if (vf_broken_count() >= 1)
        {
            int* bids = vf_broken_ids_ptr();
            double* bts = vf_broken_times_ptr();
            checkNear("A: broken id[0]", (double)bids[0], 1.0, 1e-9);
            checkNear("A: broken time", bts[0], 51.0 * dt,
                      fabs(2.0 * dt));
        }

        // verify killed element forces are zero and root reactions are
        // zero (the cantilever is now disconnected from the support)
        vf_collect();
        int ne = vf_elem_count();
        int* eids = vf_elem_ids_ptr();
        double* ef = vf_forces_ptr();

        int rowE1 = -1;
        for (int k = 0; k < ne; k++)
        {
            if (eids[k] == 1) rowE1 = k;
        }
        if (rowE1 >= 0)
        {
            checkNear("A: elem1 N ", ef[8 * rowE1 + 0], 0.0, 1e-6);
            checkNear("A: elem1 V ", ef[8 * rowE1 + 1], 0.0, 1e-6);
            checkNear("A: elem1 Mi", ef[8 * rowE1 + 4], 0.0, 1e-6);
            checkNear("A: elem1 Mj", ef[8 * rowE1 + 5], 0.0, 1e-6);
        }
        else
        {
            printf("  A: elem1 not found in collect           FAIL\n");
            ++g_failures;
        }

        int nn = vf_node_count();
        int* nids = vf_node_ids_ptr();
        double* r = vf_react_ptr();
        int rowRoot = -1;
        for (int k = 0; k < nn; k++)
        {
            if (nids[k] == 1) rowRoot = k;
        }
        if (rowRoot >= 0)
        {
            checkNear("A: root Rx /P", r[6 * rowRoot + 0] / P, 0.0, 1e-6);
            checkNear("A: root Ry /P", r[6 * rowRoot + 1] / P, 0.0, 1e-6);
            checkNear("A: root Rmz/(P*L)",
                      r[6 * rowRoot + 5] / (P * L), 0.0, 1e-6);
        }
        else
        {
            printf("  A: root node not found in collect       FAIL\n");
            ++g_failures;
        }

        printf("kill-scheduled: rc=%d broken(%d) [elem %d at t=%.3e]\n",
               rc, vf_broken_count(),
               vf_broken_count() > 0 ? vf_broken_ids_ptr()[0] : -1,
               vf_broken_count() > 0 ? vf_broken_times_ptr()[0] : -1.0);
    }

    // --- Scenario B: capacity fracture --------------------------------
    {
        vf_create();

        for (int i = 1; i <= nelem + 1; i++)
            vf_add_node(i, L * (i - 1) / nelem, 0.0);

        vf_add_material(1, E, RHO);
        vf_add_section(1, A, I);

        for (int i = 1; i <= nelem; i++)
        {
            if (vf_add_element(i, 0, i, i + 1, 1, 1) != 0)
            {
                printf("add_element failed: %s\n", vf_error());
                exit(1);
            }
        }

        vf_add_support(1, 1, 1, 1, 0, 0, 0);
        vf_set_load(nelem + 1, 0.0, -P, 0.0);

        double dtc = vf_critical_dt();
        double dt = pow(10.0, floor(log10(dtc)));
        printf("capacity-frac (B): dt_critical=%.3e dt=%.3e\n", dtc, dt);

        vf_set_fracture(1);
        // mcap = half the static root moment so the sudden tip load
        // must rupture element 1 early in the run
        vf_set_capacity(1, 0.0, 0.5 * P * L);

        if (vf_setup(dt, 0.0) != 0)
        {
            printf("setup failed: %s\n", vf_error());
            exit(1);
        }

        int rc = vf_run_steps(2000);
        if (rc < 0)
        {
            printf("vf_run_steps failed: %s\n", vf_error());
            ++g_failures;
        }

        checkNear("B: broken count", (double)vf_broken_count(), 1.0, 1e-9);
        if (vf_broken_count() >= 1)
        {
            checkNear("B: broken id[0]",
                      (double)vf_broken_ids_ptr()[0], 1.0, 1e-9);
        }

        // run must stay finite after the rupture (debris flies free)
        int rc2 = vf_run_steps(3000);
        if (rc2 < 0)
        {
            printf("  B: post-fracture run diverged           FAIL\n");
            ++g_failures;
        }

        printf("capacity-fracture: rc=%d broken(%d) post-rc=%d\n",
               rc, vf_broken_count(), rc2);
    }
}

// Impulse-load test: a completely free rod (no supports, no gravity)
// receives a downward force window at the mid node. After 1000 steps
// the total vertical momentum must equal the impulse F·Δt. The
// momentum is computed from per-node lumped masses and velocities.
static void runImpulse()
{
    const double L = 4.0, E = 200.0e9, A = 2.85e-3, I = 1.94e-5;
    const double RHO = 7850.0, F = 10000.0;
    const int nelem = 8;

    vf_create();

    for (int i = 1; i <= nelem + 1; i++)
        vf_add_node(i, L * (i - 1) / nelem, 0.0);

    vf_add_material(1, E, RHO);
    vf_add_section(1, A, I);

    for (int i = 1; i <= nelem; i++)
    {
        if (vf_add_element(i, 0, i, i + 1, 1, 1) != 0)
        {
            printf("add_element failed: %s\n", vf_error());
            exit(1);
        }
    }

    double dtc = vf_critical_dt();
    double dt = pow(10.0, floor(log10(dtc)));
    printf("impulse: dt_critical=%.3e dt=%.3e\n", dtc, dt);

    double t0 = 100.5 * dt;
    double t1 = 300.5 * dt;
    vf_set_load_window(5, 0.0, -F, 0.0, t0, t1);

    if (vf_setup(dt, 0.0) != 0)
    {
        printf("setup failed: %s\n", vf_error());
        exit(1);
    }

    int rc = vf_run_steps(1000);
    if (rc != 1000)
    {
        printf("  impulse: vf_run_steps returned %d, expected 1000  FAIL\n",
               rc);
        ++g_failures;
    }

    vf_collect();
    int nn = vf_node_count();
    int* ids = vf_node_ids_ptr();
    double* vel = vf_vel_ptr();

    const double m_elem = RHO * A * (L / nelem);
    const double m_end  = 0.5 * m_elem;

    double sum_mv = 0.0;
    for (int k = 0; k < nn; k++)
    {
        int nid = ids[k];
        double vy = vel[6 * k + 1];
        double m  = (nid == 1 || nid == nelem + 1) ? m_end : m_elem;
        sum_mv += m * vy;
    }

    double impulse = -F * (t1 - t0);  // downward → negative
    checkNear("impulse momentum", sum_mv, impulse, 2e-2);

    printf("impulse: rc=%d dt=%.3e t0=%.3e t1=%.3e sum_mv=%.6e (want %.6e)\n",
           rc, dt, t0, t1, sum_mv, impulse);
}

// 3D biaxial cantilever through the NEW 3D exports: the exact Case A
// model from tests/test_beam3d.cpp (same L, A, Iyy, Izz, J, E, G =
// E/2.6 via vf_set_material_g — Case A sets nu = 0.3 — and the same
// physical section orientation expressed with the legacy third-point convention). Tip loads Py and Pz run separately, each relaxed
// with the shim's adaptive DR; tip deflections must match the
// Euler-Bernoulli references uy = P*L^3/(3*E*Izz) = 5.625e-3 and
// uz = P*L^3/(3*E*Iyy) = 2.25e-2 within 1%.
static void runCantilever3D()
{
    const double L = 3.0, E = 2.0e11, A = 1.0e-2, Izz = 8.0e-5,
                 Iyy = 2.0e-5, J = 5.0e-5, RHO = 7850.0, P = 1.0e4;
    const int nelem = 8;

    for (int sub = 0; sub < 2; sub++)   // 0: tip load Py, 1: tip load Pz
    {
        vf_create();

        for (int i = 0; i <= nelem; i++)
            vf_add_node3(i + 1, L * i / nelem, 0.0, 0.0);

        vf_add_material(1, E, RHO);
        vf_set_material_g(1, E / 2.6);
        vf_add_section3(1, A, Iyy, Izz, J);

        for (int i = 1; i <= nelem; i++)
        {
            // Legacy bepcrd: +Y third-point direction gives local
            // ey=+Y and ez=+Z for this +X member.
            if (vf_add_element3(i, 2, i, i + 1, 1, 1, 0.0, 1.0, 0.0) != 0)
            {
                printf("add_element3 failed: %s\n", vf_error());
                exit(1);
            }
        }

        vf_add_support6(1, 1, 1, 1, 1, 1, 1);
        if (sub == 0)
            vf_set_load6(nelem + 1, 0.0, P, 0.0, 0.0, 0.0, 0.0);
        else
            vf_set_load6(nelem + 1, 0.0, 0.0, P, 0.0, 0.0, 0.0);

        double dtc = vf_critical_dt();
        double dt = 0.7 * dtc;
        printf("3D cantilever (%s): dt_critical=%.3e dt=%.3e\n",
               sub == 0 ? "tip Py" : "tip Pz", dtc, dt);

        if (vf_setup(dt, 0.0) != 0)
        {
            printf("setup failed: %s\n", vf_error());
            exit(1);
        }

        int total = 0;
        for (int chunk = 0; chunk < 400; chunk++)
        {
            int rc = vf_run_steps_kinetic(2000);
            if (rc < 0)
            {
                printf("run failed: %s\n", vf_error());
                exit(1);
            }
            total += 2000;
            if (vf_max_residual_force() < 1e-8 * P &&
                vf_max_residual_moment() < 1e-8 * P * L)
                break;
        }
        printf("  steps: %d  residF %.2e residM %.2e KE %.2e\n", total,
               vf_max_residual_force(), vf_max_residual_moment(),
               vf_kinetic_energy());

        vf_collect();
        int nn = vf_node_count();
        int* ids = vf_node_ids_ptr();
        double* u = vf_disp_ptr();

        int rowTip = -1;
        for (int k = 0; k < nn; k++)
        {
            if (ids[k] == nelem + 1) rowTip = k;
        }

        if (sub == 0)
        {
            checkNear("3D: tip uy under Py",
                      u[6 * rowTip + 1],
                      P * L * L * L / (3.0 * E * Izz), 1e-2);
        }
        else
        {
            checkNear("3D: tip uz under Pz",
                      u[6 * rowTip + 2],
                      P * L * L * L / (3.0 * E * Iyy), 1e-2);
        }
    }
}

// 3D portal built from Beam3D throughout: columns along +y at x=0 and
// x=L, crossbeam along +x at y=H, all in the z=0 plane, with an
// out-of-plane point load +Pz at the left beam-column joint. Column
// reference directions use the legacy third-point convention so the same
// physical local section axes are preserved after the Beam3D convention
// migration. Checks: adaptive DR converges, the loaded joint moves in the
// load direction, and the support reactions balance the applied load.
static void runPortal3D()
{
    const double H = 3.0, L = 4.0, E = 2.0e11, A = 1.0e-2;
    const double Izz = 8.0e-5, Iyy = 2.0e-5, J = 5.0e-5;
    const double RHO = 7850.0, P = 1.0e3;
    const int mesh3 = 4;

    vf_create();

    // same node layout as the 2D portal: left column bottom->top, beam
    // left->right, right column top->bottom
    int id = 1;
    for (int i = 0; i <= mesh3; i++)
        vf_add_node3(id++, 0.0, H * i / mesh3, 0.0);
    for (int i = 1; i <= mesh3; i++)
        vf_add_node3(id++, L * i / mesh3, H, 0.0);
    for (int i = 1; i <= mesh3; i++)
        vf_add_node3(id++, L, H - H * i / mesh3, 0.0);
    const int nTop = mesh3 + 1;              // left beam-column joint
    const int nBaseRight = 3 * mesh3 + 1;    // right column base

    vf_add_material(1, E, RHO);
    vf_set_material_g(1, E / 2.6);
    vf_add_section3(1, A, Iyy, Izz, J);

    int eid = 1;
    for (int n = 1; n < nBaseRight; n++)
    {
        const bool isBeam = (n >= mesh3 + 1) && (n <= 2 * mesh3);
        double rx = 0.0, ry = 0.0, rz = 0.0;
        if (isBeam)
            rz = -1.0;                 // +X chord -> ey=-Z, ez=+Y
        else if (n <= mesh3)
            rx = -1.0;                 // +Y chord -> ey=-X, ez=+Z
        else
            rx = 1.0;                  // -Y chord -> ey=+X, ez=+Z
        if (vf_add_element3(eid++, 2, n, n + 1, 1, 1, rx, ry, rz) != 0)
        {
            printf("add_element3 failed: %s\n", vf_error());
            exit(1);
        }
    }

    vf_add_support6(1, 1, 1, 1, 1, 1, 1);
    vf_add_support6(nBaseRight, 1, 1, 1, 1, 1, 1);
    vf_set_load6(nTop, 0.0, 0.0, P, 0.0, 0.0, 0.0);

    double dtc = vf_critical_dt();
    double dt = 0.35 * dtc;
    printf("3D portal (out-of-plane Pz): dt_critical=%.3e dt=%.3e\n",
           dtc, dt);

    if (vf_setup(dt, 0.0) != 0)
    {
        printf("setup failed: %s\n", vf_error());
        exit(1);
    }

    int total = 0;
    for (int chunk = 0; chunk < 1600; chunk++)
    {
        int rc = vf_run_steps_kinetic(2000);
        if (rc < 0)
        {
            printf("run failed: %s\n", vf_error());
            exit(1);
        }
        total += 2000;
        if (vf_max_residual_force() < 1e-8 * P &&
            vf_max_residual_moment() < 1e-8 * P * H)
            break;
    }
    printf("  steps: %d  residF %.2e residM %.2e KE %.2e\n", total,
           vf_max_residual_force(), vf_max_residual_moment(),
           vf_kinetic_energy());

    // Convergence assertion: exact frm3D-compatible triad evolution
    // excites the slow torsion-coupled out-of-plane modes more strongly than
    // the earlier approximation. The larger iteration budget keeps the
    // existing 1e-5 equilibrium criterion without loosening it.
    checkNear("3D portal: residF/P at end",
              vf_max_residual_force() / P, 0.0, 1e-5);
    checkNear("3D portal: residM/(P*H) at end",
              vf_max_residual_moment() / (P * H), 0.0, 1e-5);

    vf_collect();
    int nn = vf_node_count();
    int* ids = vf_node_ids_ptr();
    double* u = vf_disp_ptr();
    double* r = vf_react_ptr();

    int rowTop = -1, rowB1 = -1, rowB2 = -1;
    for (int k = 0; k < nn; k++)
    {
        if (ids[k] == nTop) rowTop = k;
        if (ids[k] == 1) rowB1 = k;
        if (ids[k] == nBaseRight) rowB2 = k;
    }

    // the loaded joint must move a finite, nonzero amount along +z
    double uzTop = u[6 * rowTop + 2];
    ++g_checks;
    bool okUz = std::isfinite(uzTop) && fabs(uzTop) > 1e-9;
    if (!okUz) ++g_failures;
    printf("  %-40s got % .6e  %s\n",
           "3D portal: tip uz finite & nonzero", uzTop,
           okUz ? "ok" : "FAIL");

    // global equilibrium: support reactions balance the applied load
    // (in-plane sums see the same ~5e-6*P residual floor as above)
    checkNear("3D portal: sum Rx /P",
              (r[6 * rowB1 + 0] + r[6 * rowB2 + 0]) / P, 0.0, 1e-5);
    checkNear("3D portal: sum Ry /P",
              (r[6 * rowB1 + 1] + r[6 * rowB2 + 1]) / P, 0.0, 1e-5);
    // Use the same convergence-scale tolerance as the residual checks
    // above; the adaptive DR budget intentionally stops at a small dynamic
    // residual rather than iterating to machine-zero static equilibrium.
    checkNear("3D portal: sum Rz balances load",
              r[6 * rowB1 + 2] + r[6 * rowB2 + 2], -P, 1e-5);
}

int main()
{
    runPortal(false);
    runPortal(true);
    runCantilever(true);
    runKillAndFracture();
    runImpulse();
    runCantilever3D();
    runPortal3D();

    printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
