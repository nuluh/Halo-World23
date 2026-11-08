// Modern C++ stage probe for the frozen legacy Fortran oracle fixtures.
// Emits key,value CSV; validation/fortran/golden/compare_cpp.py owns tolerance
// policy and reports the earliest divergent numerical stage.
#include <openvfife/element/frame.hpp>
#include <openvfife/material/uniaxial.hpp>
#include <openvfife/core/particle.hpp>
#include <openvfife/section/section.hpp>

#include <array>
#include <iomanip>
#include <iostream>
#include <string>

using openvfife::StdArray3d;
using openvfife::StdArray6d;
using openvfife::Dof;
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
    // A, Sy, Sz, SHy, SHz, CGy, CGz, Iyy, Izz, Iyz, J
    CustomSectionParas sp = {
        1.0e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        2.0e-5, 3.0e-5, 0.0, 1.0e-5
    };
    return sp;
}

static void emitBasis(const string& prefix, const Beam3D& beam)
{
    const Eigen::Vector3d& ex = beam.ex();
    const Eigen::Vector3d& ey = beam.ey();
    const Eigen::Vector3d& ez = beam.ez();
    const Eigen::Vector3d basis[3] = {ex, ey, ez};
    int key = 1;
    for (int axis = 0; axis < 3; ++axis)
    {
        for (int component = 0; component < 3; ++component)
        {
            emit(prefix + std::to_string(key++), basis[axis](component));
        }
    }
}

static void emitRmas(const string& prefix, const Eigen::Matrix3d& m)
{
    // Legacy rmas order: Ixx, Iyy, Izz, Iyz, Ixz, Ixy.
    emit(prefix + "1", m(0, 0));
    emit(prefix + "2", m(1, 1));
    emit(prefix + "3", m(2, 2));
    emit(prefix + "4", m(1, 2));
    emit(prefix + "5", m(0, 2));
    emit(prefix + "6", m(0, 1));
}

static void basicStages()
{
    Particle p1(1, 0.0, 0.0, 0.0);
    Particle p2(2, 2.0, 0.0, 0.0);
    LinearElastic mat(1);
    makeOracleMaterial(mat);
    CustomSectionParas sp = oracleSection();
    CustomSection section(1, sp);
    const StdArray3d ref = {0.0, 1.0, 0.0};
    Beam3D beam(1, p1, p2, mat, section, ref);

    emitBasis("basis.", beam);

    const Mass& m1 = p1.mass();
    const Mass& m2 = p2.mass();
    emit("tmass.1", m1.lumped_mass[0]);
    emit("tmass.2", m1.lumped_mass[1]);
    emit("tmass.3", m1.lumped_mass[2]);
    emit("tmass.4", m2.lumped_mass[0]);
    emit("tmass.5", m2.lumped_mass[1]);
    emit("tmass.6", m2.lumped_mass[2]);
    emitRmas("rmas_node1.", m1.rotation_mass_matrix);
    emitRmas("rmas_node2.", m2.rotation_mass_matrix);
}

static void frmCase(const string& name, const std::array<double, 12>& dd)
{
    Particle p1(1, 0.0, 0.0, 0.0);
    Particle p2(2, 2.0, 0.0, 0.0);
    LinearElastic mat(1);
    makeOracleMaterial(mat);
    CustomSectionParas sp = oracleSection();
    CustomSection section(1, sp);
    const StdArray3d ref = {0.0, 1.0, 0.0};
    Beam3D beam(1, p1, p2, mat, section, ref);

    // Legacy Dd layout is all translations first, then all rotations:
    // [ui,vi,wi, uj,vj,wj, rxi,ryi,rzi, rxj,ryj,rzj].
    StdArray6d ui = {dd[0], dd[1], dd[2], dd[6], dd[7], dd[8]};
    StdArray6d uj = {dd[3], dd[4], dd[5], dd[9], dd[10], dd[11]};
    p1.setDisplace(ui);
    p2.setDisplace(uj);
    p1.updatePosition();
    p2.updatePosition();

    beam.setParticleForce();

    // C++ Particle::force stores the residual/restoring contribution, while
    // legacy elf is the positive internal-force vector used as Fext-Fint.
    const StdArray6d& fi = p1.force();
    const StdArray6d& fj = p2.force();
    for (int i = 0; i < 6; ++i)
        emit(name + ".elf." + std::to_string(i + 1), -fi[i]);
    for (int i = 0; i < 6; ++i)
        emit(name + ".elf." + std::to_string(i + 7), -fj[i]);

    emitBasis(name + ".basis.", beam);
}

static void elementStages()
{
    std::array<double, 12> dd = {};

    dd[3] = 1.0e-3;
    frmCase("axial_extension", dd);

    dd = {};
    dd[9] = 1.0e-2;
    frmCase("torsion", dd);

    dd = {};
    dd[7] = 1.0e-2;
    frmCase("bend_y_node1", dd);

    dd = {};
    dd[8] = 1.0e-2;
    frmCase("bend_z_node1", dd);

    dd = {};
    dd[3] = 5.0e-4;
    dd[4] = 2.0e-4;
    dd[5] = -3.0e-4;
    dd[9] = 4.0e-3;
    dd[10] = -6.0e-3;
    dd[11] = 8.0e-3;
    frmCase("combined", dd);

    const double theta = 1.0e-2;
    dd = {};
    dd[3] = 2.0 * (std::cos(theta) - 1.0);
    dd[4] = 2.0 * std::sin(theta);
    dd[8] = theta;
    dd[11] = theta;
    frmCase("rigid_z", dd);
}

static void emitVector(const string& prefix, const StdArray6d& values)
{
    for (int i = 0; i < 6; ++i)
        emit(prefix + std::to_string(i + 1), values[i]);
}

static void integrationStage()
{
    Particle p(1, 0.0, 0.0, 0.0);
    p.activateDof(Dof::Ux);
    p.activateDof(Dof::Uy);
    p.activateDof(Dof::Uz);
    p.activateDof(Dof::Rotx);
    p.activateDof(Dof::Roty);
    p.activateDof(Dof::Rotz);

    const StdArray6d lumped = {2.0, 4.0, 5.0, 0.0, 0.0, 0.0};
    p.setLumpedMass(lumped);

    Eigen::Matrix3d inertia;
    inertia << 2.0, 0.3, 0.1,
               0.3, 3.0, 0.2,
               0.1, 0.2, 4.0;
    p.setInertiaMass(inertia);

    const StdArray6d d0 = {0.1, -0.05, 0.02, 0.01, -0.02, 0.03};
    const StdArray6d v0 = {0.2, -0.1, 0.05, 0.02, 0.01, -0.03};
    const StdArray6d force = {10.0, -8.0, 4.0, 6.0, -5.0, 3.0};
    p.setDisplace(d0);
    p.setVelocity(v0);
    p.setForce(force);
    p.updatePosition();

    const double h = 0.1;
    const double damping = 0.4;
    p.solveFirstStep(h, damping);
    const StdArray6d dprev = p.previous_displace();
    emitVector("dprev_isolated.", dprev);

    p.solve(h, damping);
    const StdArray6d dnext = p.displace();
    emitVector("dnext_isolated.", dnext);
    emitVector("velocity_isolated.", p.velocity());
    emitVector("accel_isolated.", p.accelerate());

    StdArray6d increment = {};
    for (int i = 0; i < 6; ++i)
        increment[i] = dnext[i] - d0[i];
    emitVector("increment_isolated.", increment);
}

int main()
{
    basicStages();
    elementStages();
    integrationStage();
    return 0;
}
