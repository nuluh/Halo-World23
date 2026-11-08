#include <openvfife/element/frame.hpp>
#include <openvfife/element/experimental.hpp>

#include <cstdio>
#include <stdexcept>
#include <type_traits>

using openvfife::StdArray3d;
using openvfife::Particle;
using openvfife::BaseMaterial;
using openvfife::LinearElastic;
using openvfife::BaseSection;
using openvfife::CustomSectionParas;
using openvfife::CustomSection;
using openvfife::Link2D;
using openvfife::Link3D;
using openvfife::Beam2D;
using openvfife::Beam3D;

namespace {
int checks = 0;
int failures = 0;

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
    std::printf("  %-42s %s\n", label, ok ? "ok" : "FAIL");
}

void check_true(const char* label, bool value)
{
    ++checks;
    if (!value) ++failures;
    std::printf("  %-42s %s\n", label, value ? "ok" : "FAIL");
}
} // namespace

int main()
{
    LinearElastic material(1);
    material.setE(200.0e9);
    material.setDensity(7850.0);

    CustomSectionParas params{};
    params.A = 1.0e-2;
    params.Iyy = 2.0e-5;
    params.Izz = 3.0e-5;
    CustomSection section(1, params);

    Particle p1(1, 0.0, 0.0, 0.0);
    Particle p2(2, 2.0, 0.0, 0.0);
    Particle p_same_id(1, 2.0, 0.0, 0.0);
    Particle p_same_position(3, 0.0, 0.0, 0.0);

    static_assert(std::is_constructible_v<Link2D, int, Particle&, Particle&,
                                          const BaseMaterial&, const BaseSection&>);
    static_assert(!std::is_constructible_v<Link2D, int, Particle*, Particle*,
                                           const BaseMaterial*, const BaseSection*>);
    static_assert(std::is_constructible_v<Link3D, int, Particle&, Particle&,
                                          const BaseMaterial&, const BaseSection&>);
    static_assert(!std::is_constructible_v<Link3D, int, Particle*, Particle*,
                                           const BaseMaterial*, const BaseSection*>);
    static_assert(std::is_constructible_v<Beam2D, int, Particle&, Particle&,
                                          const BaseMaterial&, const BaseSection&>);
    static_assert(!std::is_constructible_v<Beam2D, int, Particle*, Particle*,
                                           const BaseMaterial*, const BaseSection*>);
    static_assert(std::is_constructible_v<Beam3D, int, Particle&, Particle&,
                                          const BaseMaterial&, const BaseSection&>);
    static_assert(std::is_constructible_v<Beam3D, int, Particle&, Particle&,
                                          const BaseMaterial&, const BaseSection&, StdArray3d>);
    static_assert(!std::is_constructible_v<Beam3D, int, Particle*, Particle*,
                                           const BaseMaterial*, const BaseSection*>);

    check_throws<std::invalid_argument>("reject same particle twice", [&] {
        Link2D invalid(1, p1, p1, material, section);
    });
    check_throws<std::invalid_argument>("reject duplicate particle ids", [&] {
        Link2D invalid(1, p1, p_same_id, material, section);
    });
    check_throws<std::invalid_argument>("reject zero element length", [&] {
        Link2D invalid(1, p1, p_same_position, material, section);
    });

    Link2D valid(2, p1, p2, material, section);
    check_true("valid element has two particles", valid.particles().size() == 2);
    check_true("valid element has positive length", valid.length() > 0.0);

    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
