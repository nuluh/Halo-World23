#include <openvfife/material/uniaxial.hpp>
#include <openvfife/material/experimental.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>

using openvfife::UniaxialMaterialState;
using openvfife::LinearElastic;
using openvfife::UniIdeal;
using openvfife::UniBilinear;
using openvfife::UniPower;

namespace {
int checks = 0;
int failures = 0;

void check_near(const char* label, double actual, double expected,
                double tolerance = 1.0e-12)
{
    ++checks;
    const bool ok = std::abs(actual - expected) <= tolerance;
    if (!ok) {
        ++failures;
    }
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
    if (!ok) {
        ++failures;
    }
    std::printf("  %-38s %s\n", label, ok ? "ok" : "FAIL");
}

void check_state(const UniaxialMaterialState& state,
                 double stress, double strain, double plastic_strain,
                 double back_stress, double hardening)
{
    check_near("state stress", state.stress, stress);
    check_near("state strain", state.strain, strain);
    check_near("state plastic strain", state.plastic_strain, plastic_strain);
    check_near("state back stress", state.back_stress, back_stress);
    check_near("state hardening", state.hardening, hardening);
}
} // namespace

int main()
{
    LinearElastic elastic(1);
    elastic.setE(200.0e9);
    UniaxialMaterialState elastic_state;
    const double elastic_tangent = elastic.update(1.0e-3, elastic_state);
    check_near("elastic tangent", elastic_tangent, 200.0e9);
    check_state(elastic_state, 200.0e6, 1.0e-3, 0.0, 0.0, 0.0);

    UniBilinear bilinear(2, 200.0, 20.0, 100.0, 0.5);
    UniaxialMaterialState state;

    check_near("bilinear elastic tangent", bilinear.update(0.2, state), 200.0);
    check_state(state, 40.0, 0.2, 0.0, 0.0, 0.0);

    check_near("bilinear plastic tangent", bilinear.update(0.4, state), 20.0);
    check_state(state, 102.0, 0.6, 0.09, 2.0, 2.0);

    check_near("bilinear unload tangent", bilinear.update(-0.1, state), 200.0);
    check_state(state, 82.0, 0.5, 0.09, 2.0, 2.0);

    elastic.setE(150.0);
    check_throws<std::invalid_argument>("reject negative E", [&] { elastic.setE(-1.0); });
    check_near("E unchanged after rejection", elastic.E(), 150.0);

    elastic.setG(60.0);
    check_throws<std::invalid_argument>("reject negative G", [&] { elastic.setG(-1.0); });
    check_near("G unchanged after rejection", elastic.G(), 60.0);

    elastic.setNu(0.25);
    const double accepted_g = elastic.G();
    check_throws<std::invalid_argument>("reject negative nu", [&] { elastic.setNu(-0.1); });
    check_near("nu unchanged after rejection", elastic.nu(), 0.25);
    check_near("derived G unchanged", elastic.G(), accepted_g);

    elastic.setDensity(7850.0);
    check_throws<std::invalid_argument>("reject negative density", [&] { elastic.setDensity(-1.0); });
    check_near("density unchanged after rejection", elastic.density(), 7850.0);

    elastic.setLimitedTensileStrain(0.02);
    check_throws<std::invalid_argument>("reject nonpositive tensile limit", [&] {
        elastic.setLimitedTensileStrain(0.0);
    });
    check_near("tensile limit unchanged", elastic.limited_tensile_strain(), 0.02);

    bilinear.setYieldStress(120.0);
    check_throws<std::invalid_argument>("reject nonpositive yield stress", [&] {
        bilinear.setYieldStress(0.0);
    });
    check_near("yield stress unchanged", bilinear.yield_stress(), 120.0);

    check_throws<std::invalid_argument>("UniIdeal rejects zero E", [] {
        UniIdeal invalid(10, 0.0, 100.0);
    });
    check_throws<std::invalid_argument>("UniIdeal rejects zero yield", [] {
        UniIdeal invalid(10, 200.0, 0.0);
    });
    check_throws<std::invalid_argument>("UniBilinear rejects zero E", [] {
        UniBilinear invalid(11, 0.0, 0.0, 100.0, 0.5);
    });
    check_throws<std::invalid_argument>("UniBilinear rejects negative Et", [] {
        UniBilinear invalid(11, 200.0, -1.0, 100.0, 0.5);
    });
    check_throws<std::invalid_argument>("UniBilinear rejects Et above E", [] {
        UniBilinear invalid(11, 200.0, 201.0, 100.0, 0.5);
    });
    check_throws<std::invalid_argument>("UniBilinear rejects zero yield", [] {
        UniBilinear invalid(11, 200.0, 20.0, 0.0, 0.5);
    });
    check_throws<std::invalid_argument>("UniBilinear rejects hardening mix", [] {
        UniBilinear invalid(11, 200.0, 20.0, 100.0, 1.1);
    });
    check_throws<std::invalid_argument>("UniPower rejects zero E", [] {
        UniPower invalid(12, 0.0, 0.2, 100.0, 0.5);
    });
    check_throws<std::invalid_argument>("UniPower rejects zero yield", [] {
        UniPower invalid(12, 200.0, 0.2, 0.0, 0.5);
    });
    check_throws<std::invalid_argument>("UniPower rejects hardening mix", [] {
        UniPower invalid(12, 200.0, 0.2, 100.0, -0.1);
    });

    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
