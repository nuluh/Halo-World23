#include <openvfife/numerics/quadrature.hpp>

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

using openvfife::GaussLegendre;
using openvfife::GaussLobatto;

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
    std::printf("  %-34s got % .12e want % .12e %s\n",
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
    std::printf("  %-34s %s\n", label, ok ? "ok" : "FAIL");
}
} // namespace

int main()
{
    GaussLegendre q2(2);
    const std::vector<double> one_d{1.0, 1.0};
    check_near("1D constant integral", q2.integrate1d(one_d, -1.0, 1.0), 2.0);

    const std::vector<std::vector<double>> two_d{
        {1.0, 1.0},
        {1.0, 1.0},
    };
    check_near("2D constant integral",
               q2.integrate2d(two_d, -1.0, 1.0, -1.0, 1.0), 4.0);

    const std::vector<std::vector<std::vector<double>>> three_d(
        2, std::vector<std::vector<double>>(2, std::vector<double>(2, 1.0)));
    check_near("3D constant integral",
               q2.integrate3d(three_d, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0), 8.0);

    check_throws<std::invalid_argument>("Legendre rejects zero order", [] {
        GaussLegendre invalid(0);
    });
    check_throws<std::invalid_argument>("Legendre rejects negative order", [] {
        GaussLegendre invalid(-1);
    });

    check_throws<std::out_of_range>("x rejects invalid index", [&] { q2.x(2); });
    check_throws<std::out_of_range>("w rejects negative index", [&] { q2.w(-1); });
    check_throws<std::invalid_argument>("1D rejects short input", [&] {
        q2.integrate1d(std::vector<double>{1.0}, -1.0, 1.0);
    });
    check_throws<std::invalid_argument>("2D rejects empty input", [&] {
        q2.integrate2d({}, -1.0, 1.0, -1.0, 1.0);
    });
    check_throws<std::invalid_argument>("3D rejects empty plane", [&] {
        q2.integrate3d({{}}, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    });
    check_throws<std::invalid_argument>("Legendre rejects order above 10", [] {
        GaussLegendre invalid(11);
    });
    check_throws<std::invalid_argument>("Lobatto rejects order below 2", [] {
        GaussLobatto invalid(1);
    });
    check_throws<std::invalid_argument>("Lobatto rejects order above 10", [] {
        GaussLobatto invalid(11);
    });

    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
