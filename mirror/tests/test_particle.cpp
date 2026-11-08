#include <openvfife/core/particle.hpp>

#include <cmath>
#include <cstdio>
#include <stdexcept>

using openvfife::StdArray6d;
using openvfife::Dof;
using openvfife::Particle;

namespace {
int checks = 0;
int failures = 0;

void check_near(const char* label, double actual, double expected,
                double tolerance = 0.0)
{
    ++checks;
    const bool ok = std::abs(actual - expected) <= tolerance;
    if (!ok) ++failures;
    std::printf("  %-38s got % .12e want % .12e %s\n",
                label, actual, expected, ok ? "ok" : "FAIL");
}

template <typename Exception, typename F>
void check_throws(const char* label, F&& fn)
{
    ++checks;
    bool ok = false;
    try {
        fn();
    } catch (const Exception&) {
        ok = true;
    } catch (...) {
    }
    if (!ok) ++failures;
    std::printf("  %-38s %s\n", label, ok ? "ok" : "FAIL");
}

void check_zero_state(const char* prefix, const Particle& p)
{
    for (std::size_t i = 0; i < 6; ++i) {
        check_near((std::string(prefix) + " coordinate").c_str(), p.coordinate()[i], 0.0);
        check_near((std::string(prefix) + " force").c_str(), p.force()[i], 0.0);
        check_near((std::string(prefix) + " displacement").c_str(), p.displace()[i], 0.0);
        check_near((std::string(prefix) + " velocity").c_str(), p.velocity()[i], 0.0);
        check_near((std::string(prefix) + " acceleration").c_str(), p.accelerate()[i], 0.0);
        check_near((std::string(prefix) + " current position").c_str(), p.current_position()[i], 0.0);
        check_near((std::string(prefix) + " previous position").c_str(), p.previous_position()[i], 0.0);
    }
}
} // namespace

int main()
{
    Particle origin(1);
    check_zero_state("origin", origin);
    check_near("origin translate mass norm", origin.mass().translate_mass_matrix.norm(), 0.0);
    check_near("origin rotation mass norm", origin.mass().rotation_mass_matrix.norm(), 0.0);

    Particle typed_dof(7);
    typed_dof.activateDof(Dof::Ux);
    check_near("typed dof activates Ux", typed_dof.dof_index()[0] ? 1.0 : 0.0, 1.0);
    typed_dof.deactivateDof(Dof::Ux);
    check_near("typed dof deactivates Ux", typed_dof.dof_index()[0] ? 1.0 : 0.0, 0.0);

    check_throws<std::invalid_argument>("reject invalid typed dof", [&] {
        typed_dof.activateDof(static_cast<Dof>(255));
    });

    Particle planar(2, 1.25, -3.5);
    check_near("planar x", planar.coordinate()[0], 1.25);
    check_near("planar y", planar.coordinate()[1], -3.5);
    check_near("planar z", planar.coordinate()[2], 0.0);
    check_near("planar current x", planar.current_position()[0], 1.25);
    check_near("planar previous y", planar.previous_position()[1], -3.5);

    Particle spatial(3, 1.0, 2.0, 3.0);
    check_near("spatial x", spatial.coordinate()[0], 1.0);
    check_near("spatial y", spatial.coordinate()[1], 2.0);
    check_near("spatial z", spatial.coordinate()[2], 3.0);
    check_near("spatial current z", spatial.current_position()[2], 3.0);
    check_near("spatial previous z", spatial.previous_position()[2], 3.0);

    StdArray6d accepted_mass{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    spatial.setLumpedMass(accepted_mass);
    check_near("accepted lumped mass", spatial.mass().lumped_mass[2], 3.0);

    StdArray6d invalid_mass = accepted_mass;
    invalid_mass[2] = -1.0;
    check_throws<std::invalid_argument>("reject negative lumped mass", [&] {
        spatial.setLumpedMass(invalid_mass);
    });
    check_near("set rejection preserves mass", spatial.mass().lumped_mass[2], 3.0);

    StdArray6d invalid_delta{};
    invalid_delta[1] = -3.0;
    check_throws<std::invalid_argument>("reject negative resulting mass", [&] {
        spatial.addLumpedMass(invalid_delta);
    });
    check_near("add rejection preserves mass", spatial.mass().lumped_mass[1], 2.0);

    StdArray6d safe_delta{};
    safe_delta[1] = -1.0;
    spatial.addLumpedMass(safe_delta);
    check_near("allow safe mass reduction", spatial.mass().lumped_mass[1], 1.0);

    check_throws<std::invalid_argument>("reject zero first-step dt", [&] {
        spatial.solveFirstStep(0.0, 0.2);
    });
    check_throws<std::invalid_argument>("reject negative solve dt", [&] {
        spatial.solve(-1.0e-3, 0.2);
    });

    check_throws<std::invalid_argument>("reject negative damping", [&] {
        spatial.solve(1.0e-3, -0.1);
    });

    Particle translate_matrix(5);
    translate_matrix.activateDof(Dof::Ux);
    StdArray6d one_mass{};
    one_mass[0] = 1.0;
    translate_matrix.setLumpedMass(one_mass);
    translate_matrix.setTranslateMass(Eigen::Matrix3d::Identity());
    check_throws<std::logic_error>("reject unsupported translation mass matrix", [&] {
        translate_matrix.solveFirstStep(1.0e-3, 0.2);
    });

    Particle singular_rotation(6);
    singular_rotation.setInertiaMass(Eigen::Matrix3d::Zero());
    check_throws<std::domain_error>("reject singular rotation mass matrix", [&] {
        singular_rotation.solveFirstStep(1.0e-3, 0.2);
    });

    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
