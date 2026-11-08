#include <openvfife/section/section.hpp>
#include <openvfife/section/experimental.hpp>

#include <cmath>
#include <cstdio>
#include <stdexcept>

using openvfife::CustomSectionParas;
using openvfife::CustomSection;
using openvfife::AngleSteel;

namespace {
int checks = 0;
int failures = 0;

void check_true(const char* label, bool value)
{
    ++checks;
    if (!value) ++failures;
    std::printf("  %-38s %s\n", label, value ? "ok" : "FAIL");
}

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
} // namespace

int main()
{
    CustomSectionParas polar{};
    polar.A = 2.0;
    polar.Iyy = 3.0;
    polar.Izz = 4.0;
    CustomSection fallback(1, polar);
    check_true("custom type label", fallback.type() == "CustomSection");
    check_near("custom area", fallback.A(), 2.0);
    check_near("polar J fallback", fallback.J(), 7.0);

    auto explicit_j = polar;
    explicit_j.J = 5.5;
    CustomSection explicit_section(2, explicit_j);
    check_near("explicit J", explicit_section.J(), 5.5);

    auto negative = polar;
    negative.A = -1.0;
    check_throws<std::invalid_argument>("reject negative custom area", [&] {
        CustomSection invalid(3, negative);
    });

    AngleSteel angle(4, 1.0, 2.0, 0.2, 0.1);
    check_true("valid angle area positive", angle.A() > 0.0);
    check_throws<std::invalid_argument>("reject angle t1 above w2", [] {
        AngleSteel invalid(5, 1.0, 0.1, 0.2, 0.1);
    });
    check_throws<std::invalid_argument>("reject angle t2 above w1", [] {
        AngleSteel invalid(6, 0.1, 1.0, 0.1, 0.2);
    });

    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
