// Multi-element C++ progression probe for legacy 3D-frame parity.
// This intentionally mirrors frame_progression_oracle.f90 case topology and
// output ordering, while using only public OpenVFIFE APIs.
#include <openvfife/element/frame.hpp>
#include <openvfife/material/uniaxial.hpp>
#include <openvfife/core/particle.hpp>
#include <openvfife/section/section.hpp>

#include <array>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using openvfife::StdArray3d;
using openvfife::StdArray6d;
using openvfife::DOF;
using openvfife::Mass;
using openvfife::Particle;
using openvfife::LinearElastic;
using openvfife::CustomSectionParas;
using openvfife::CustomSection;
using openvfife::Beam3D;

using std::cout;
using std::setprecision;
using std::scientific;
using std::string;
using std::unique_ptr;
using std::vector;

static void emit(const string& key, double value)
{
    cout << key << "," << scientific << setprecision(17) << value << "\n";
}

static void makeOracleMaterial(LinearElastic& mat)
{
    mat.setE(2.0e11);
    mat.setG(8.0e10);
    mat.setDensity(7850.0);
}

static CustomSectionParas oracleSection()
{
    CustomSectionParas sp = {
        1.0e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        2.0e-5, 3.0e-5, 0.0, 1.0e-5
    };
    return sp;
}

static void emitBasis(const string& prefix, const Beam3D& beam)
{
    const Eigen::Vector3d basis[3] = {beam.ex(), beam.ey(), beam.ez()};
    int key = 1;
    for (int axis = 0; axis < 3; ++axis)
        for (int component = 0; component < 3; ++component)
            emit(prefix + std::to_string(key++), basis[axis](component));
}

static void emitRmas(const string& prefix, const Eigen::Matrix3d& m)
{
    emit(prefix + "1", m(0, 0));
    emit(prefix + "2", m(1, 1));
    emit(prefix + "3", m(2, 2));
    emit(prefix + "4", m(1, 2));
    emit(prefix + "5", m(0, 2));
    emit(prefix + "6", m(0, 1));
}

struct FrameCase
{
    string name;
    vector<StdArray3d> coords;
    vector<std::array<int, 2> > connectivity;  // zero-based node indices
    vector<StdArray3d> refs;
    vector<double> legacy_dd;                  // translations first, rotations second
};

static void emitStaticFrame(const FrameCase& c)
{
    LinearElastic mat(1);
    makeOracleMaterial(mat);
    CustomSectionParas sp = oracleSection();
    CustomSection section(1, sp);

    vector<unique_ptr<Particle> > particles;
    for (std::size_t i = 0; i < c.coords.size(); ++i)
    {
        particles.push_back(std::make_unique<Particle>(
            static_cast<int>(i + 1),
            c.coords[i][0], c.coords[i][1], c.coords[i][2]));
    }

    vector<unique_ptr<Beam3D> > beams;
    for (std::size_t e = 0; e < c.connectivity.size(); ++e)
    {
        beams.push_back(std::make_unique<Beam3D>(
            static_cast<int>(e + 1),
            *particles[c.connectivity[e][0]],
            *particles[c.connectivity[e][1]],
            mat, section, c.refs[e]));
    }

    const std::size_t n = particles.size();
    const std::size_t meqt = 3 * n;
    if (c.legacy_dd.size() != 6 * n)
        throw std::runtime_error("legacy_dd size does not match frame node count");

    for (std::size_t node = 0; node < n; ++node)
    {
        StdArray6d d = {
            c.legacy_dd[3 * node + 0],
            c.legacy_dd[3 * node + 1],
            c.legacy_dd[3 * node + 2],
            c.legacy_dd[meqt + 3 * node + 0],
            c.legacy_dd[meqt + 3 * node + 1],
            c.legacy_dd[meqt + 3 * node + 2]
        };
        particles[node]->setDisplace(d);
        particles[node]->updatePosition();
        particles[node]->clearForce();
    }

    for (std::size_t e = 0; e < beams.size(); ++e)
        beams[e]->setParticleForce();

    for (std::size_t e = 0; e < beams.size(); ++e)
        emitBasis(c.name + ".element" + std::to_string(e + 1) + ".basis.", *beams[e]);

    // Legacy global ordering: all translation equations by node first, then rotations.
    int eq = 1;
    for (std::size_t node = 0; node < n; ++node)
    {
        const StdArray6d& f = particles[node]->force();
        for (int k = 0; k < 3; ++k)
            emit(c.name + ".fint." + std::to_string(eq++), -f[k]);
    }
    for (std::size_t node = 0; node < n; ++node)
    {
        const StdArray6d& f = particles[node]->force();
        for (int k = 3; k < 6; ++k)
            emit(c.name + ".fint." + std::to_string(eq++), -f[k]);
    }

    eq = 1;
    for (std::size_t node = 0; node < n; ++node)
    {
        const Mass& m = particles[node]->mass();
        for (int k = 0; k < 3; ++k)
            emit(c.name + ".tmass." + std::to_string(eq++), m.lumped_mass[k]);
        emitRmas(
            c.name + ".rmas_node" + std::to_string(node + 1) + ".",
            m.rotation_mass_matrix);
    }
}

static FrameCase continuousBeam()
{
    FrameCase c;
    c.name = "continuous_beam";
    c.coords = {
        StdArray3d{{0.0, 0.0, 0.0}},
        StdArray3d{{2.0, 0.0, 0.0}},
        StdArray3d{{4.0, 0.0, 0.0}}
    };
    c.connectivity = {
        std::array<int,2>{{0,1}},
        std::array<int,2>{{1,2}}
    };
    c.refs = {
        StdArray3d{{0.0, 1.0, 0.0}},
        StdArray3d{{0.0, 1.0, 0.0}}
    };
    c.legacy_dd.assign(18, 0.0);
    c.legacy_dd[5] = -5.0e-4;   // legacy dd(6)
    c.legacy_dd[7] =  1.0e-3;   // legacy dd(8)
    c.legacy_dd[13] = 1.0e-3;   // legacy dd(14)
    c.legacy_dd[17] = 2.0e-3;   // legacy dd(18)
    return c;
}

static FrameCase planarPortal()
{
    FrameCase c;
    c.name = "planar_portal";
    c.coords = {
        StdArray3d{{0.0, 0.0, 0.0}},
        StdArray3d{{0.0, 0.0, 3.0}},
        StdArray3d{{4.0, 0.0, 3.0}},
        StdArray3d{{4.0, 0.0, 0.0}}
    };
    c.connectivity = {
        std::array<int,2>{{0,1}},
        std::array<int,2>{{1,2}},
        std::array<int,2>{{2,3}}
    };
    // The Fortran oracle uses third points at (0,1,0), (0,1,0), (4,1,0).
    // Beam3D accepts the same node-1-relative direction vectors here.
    c.refs = {
        StdArray3d{{0.0, 1.0, 0.0}},
        StdArray3d{{0.0, 1.0, 0.0}},
        StdArray3d{{0.0, 1.0, 0.0}}
    };
    c.legacy_dd.assign(24, 0.0);
    c.legacy_dd[3] = 1.0e-3;    // legacy dd(4)
    c.legacy_dd[6] = 1.2e-3;    // legacy dd(7)
    c.legacy_dd[8] = 2.0e-4;    // legacy dd(9): keeps every portal element on a defined chord-rotation path
    c.legacy_dd[16] = 2.0e-3;   // legacy dd(17)
    c.legacy_dd[19] = -1.5e-3;  // legacy dd(20)
    return c;
}

static FrameCase spatialFrame()
{
    FrameCase c;
    c.name = "spatial_frame";
    c.coords = {
        StdArray3d{{0.0, 0.0, 0.0}},
        StdArray3d{{2.0, 0.0, 0.0}},
        StdArray3d{{2.0, 2.0, 0.0}},
        StdArray3d{{2.0, 2.0, 2.0}}
    };
    c.connectivity = {
        std::array<int,2>{{0,1}},
        std::array<int,2>{{1,2}},
        std::array<int,2>{{2,3}}
    };
    c.refs = {
        StdArray3d{{0.0, 1.0, 0.0}},
        StdArray3d{{0.0, 0.0, 1.0}},
        StdArray3d{{1.0, 0.0, 0.0}}
    };
    c.legacy_dd.assign(24, 0.0);
    c.legacy_dd[6] = 1.0e-4;
    c.legacy_dd[7] = -2.0e-4;
    c.legacy_dd[8] = 1.5e-4;
    c.legacy_dd[9] = 3.0e-4;
    c.legacy_dd[10] = -4.0e-4;
    c.legacy_dd[11] = 5.0e-4;
    c.legacy_dd[18] = 1.0e-3;
    c.legacy_dd[19] = -1.5e-3;
    c.legacy_dd[20] = 2.0e-3;
    c.legacy_dd[21] = 2.0e-3;
    c.legacy_dd[22] = -3.0e-3;
    c.legacy_dd[23] = 4.0e-3;
    return c;
}

static void emitLegacyVector(
    const string& prefix,
    const vector<unique_ptr<Particle> >& particles,
    const StdArray6d& (Particle::*getter)() const noexcept)
{
    int eq = 1;
    for (std::size_t node = 0; node < particles.size(); ++node)
    {
        const StdArray6d& v = (particles[node].get()->*getter)();
        for (int k = 0; k < 3; ++k)
            emit(prefix + std::to_string(eq++), v[k]);
    }
    for (std::size_t node = 0; node < particles.size(); ++node)
    {
        const StdArray6d& v = (particles[node].get()->*getter)();
        for (int k = 3; k < 6; ++k)
            emit(prefix + std::to_string(eq++), v[k]);
    }
}

static void shortTransient()
{
    LinearElastic mat(1);
    makeOracleMaterial(mat);
    CustomSectionParas sp = oracleSection();
    CustomSection section(1, sp);

    vector<unique_ptr<Particle> > p;
    p.push_back(std::make_unique<Particle>(1, 0.0, 0.0, 0.0));
    p.push_back(std::make_unique<Particle>(2, 2.0, 0.0, 0.0));
    p.push_back(std::make_unique<Particle>(3, 4.0, 0.0, 0.0));

    const StdArray3d ref = {{0.0, 1.0, 0.0}};
    vector<unique_ptr<Beam3D> > b;
    b.push_back(std::make_unique<Beam3D>(1, *p[0], *p[1], mat, section, ref));
    b.push_back(std::make_unique<Beam3D>(2, *p[1], *p[2], mat, section, ref));

    DOF fixed;
    fixed.key = {{true, true, true, true, true, true}};
    fixed.val = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    p[0]->constraintDof(fixed);

    const double dt = 1.0e-5;
    const double damping = 0.0;
    const StdArray6d load = {{0.0, -100.0, 0.0, 0.0, 0.0, 0.0}};

    for (std::size_t i = 0; i < p.size(); ++i)
        p[i]->clearForce();
    p[1]->addForce(load);
    p[2]->addForce(load);

    // VFIFE builds d_-1 before entering the explicit time loop.
    for (std::size_t i = 0; i < p.size(); ++i)
        p[i]->solveFirstStep(dt, damping);
    p[0]->constraintDof(fixed);

    for (int step = 0; step < 4; ++step)
    {
        // Matches VFIFE 8888: solve from the force state assembled at the
        // end of the previous step, then refresh Fext and Fint.
        for (std::size_t i = 0; i < p.size(); ++i)
            p[i]->solve(dt, damping);
        p[0]->constraintDof(fixed);

        for (std::size_t i = 0; i < p.size(); ++i)
            p[i]->clearForce();
        p[1]->addForce(load);
        p[2]->addForce(load);
        for (std::size_t e = 0; e < b.size(); ++e)
            b[e]->setParticleForce();
    }

    emitLegacyVector("short_transient.d.", p, &Particle::displace);
    emitLegacyVector("short_transient.v.", p, &Particle::velocity);
    emitLegacyVector("short_transient.a.", p, &Particle::accelerate);

    // Reassemble at zero increment without Fext so the emitted value is the
    // legacy Fint vector rather than the C++ residual Fext-Fint.
    for (std::size_t i = 0; i < p.size(); ++i)
        p[i]->clearForce();
    for (std::size_t e = 0; e < b.size(); ++e)
        b[e]->setParticleForce();

    int eq = 1;
    for (std::size_t node = 0; node < p.size(); ++node)
    {
        const StdArray6d& f = p[node]->force();
        for (int k = 0; k < 3; ++k)
            emit("short_transient.fint." + std::to_string(eq++), -f[k]);
    }
    for (std::size_t node = 0; node < p.size(); ++node)
    {
        const StdArray6d& f = p[node]->force();
        for (int k = 3; k < 6; ++k)
            emit("short_transient.fint." + std::to_string(eq++), -f[k]);
    }
}

int main()
{
    emitStaticFrame(continuousBeam());
    emitStaticFrame(planarPortal());
    emitStaticFrame(spatialFrame());
    shortTransient();
    return 0;
}
