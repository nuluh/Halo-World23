#include <openvfife/solver/struct_system.hpp>
#include <openvfife/element/frame.hpp>
#include "../src/solver/filesystem_helpers.hpp"

#include <cstdio>

#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <type_traits>

using openvfife::DOF;
using openvfife::Particle;
using openvfife::BaseMaterial;
using openvfife::LinearElastic;
using openvfife::BaseSection;
using openvfife::BaseElement;
using openvfife::StructSystem;
using openvfife::makedir;
using openvfife::isDirExists;
using openvfife::isFileExists;

namespace {
int checks = 0;
int failures = 0;

void check_equal(const char* label, int actual, int expected)
{
    ++checks;
    const bool ok = actual == expected;
    if (!ok) {
        ++failures;
    }
    std::printf("  %-34s got %d want %d %s\n",
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
    std::printf("  %-34s %s\n", label, ok ? "ok" : "FAIL");
}
struct CountingMaterial : LinearElastic
{
    static inline int destroyed = 0;

    explicit CountingMaterial(int id) : LinearElastic(id) {}
    ~CountingMaterial() override { ++destroyed; }
};

} // namespace

static_assert(!std::is_copy_constructible_v<StructSystem>);
static_assert(!std::is_copy_assignable_v<StructSystem>);
static_assert(!std::is_move_constructible_v<StructSystem>);
static_assert(!std::is_move_assignable_v<StructSystem>);

int main()
{
    StructSystem system(1);
    check_equal("workdir uses current path",
                system.workdir() == std::filesystem::current_path().string(), 1);
    check_throws<std::invalid_argument>("reject zero stored dt", [&] {
        system.setTimeStep(0.0);
    });
    check_throws<std::invalid_argument>("reject negative explicit dt", [&] {
        system.solve(-1.0e-3, 0.2, false);
    });
    check_throws<std::invalid_argument>("reject negative damping coefficient", [&] {
        system.setDampCoeff(-0.1, 0.0);
    });
    check_throws<std::invalid_argument>("reject negative damping beta", [&] {
        system.setDampCoeff(0.1, -0.1);
    });
    system.addParticle(std::make_unique<Particle>(1, 0.0, 0.0, 0.0));
    check_throws<std::invalid_argument>("reject null particle owner", [&] {
        system.addParticle(std::unique_ptr<Particle>{});
    });
    check_throws<std::invalid_argument>("reject duplicate particle id", [&] {
        system.addParticle(std::make_unique<Particle>(1, 1.0, 0.0, 0.0));
    });
    check_equal("duplicate particle preserves original",
                static_cast<int>(system.particles().size()), 1);
    check_throws<std::invalid_argument>("reject null element owner", [&] {
        system.addElement(std::unique_ptr<BaseElement>{});
    });
    check_throws<std::invalid_argument>("reject null material owner", [&] {
        system.addMaterial(std::unique_ptr<BaseMaterial>{});
    });
    check_throws<std::invalid_argument>("reject null section owner", [&] {
        system.addSection(std::unique_ptr<BaseSection>{});
    });

    const DOF fixed{{true, true, false, false, false, true},
                    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    system.addConstraint(1, fixed);
    check_throws<std::invalid_argument>("reject duplicate constraint", [&] {
        system.addConstraint(1, fixed);
    });
    check_throws<std::out_of_range>("reject constraint on missing particle", [&] {
        system.addConstraint(99, fixed);
    });

    DOF changed{{true, false, false, false, false, false},
                {1.25, 0.0, 0.0, 0.0, 0.0, 0.0}};
    system.changeConstraint(1, changed);
    check_equal("change existing constraint",
                system.constraints().at(1).key[0] ? 1 : 0, 1);
    check_throws<std::out_of_range>("reject change on missing particle", [&] {
        system.changeConstraint(99, changed);
    });

    check_equal("delete existing constraint", system.deleteConstraint(1) ? 1 : 0, 1);
    check_equal("delete missing constraint", system.deleteConstraint(1) ? 1 : 0, 0);
    check_throws<std::out_of_range>("reject change without constraint", [&] {
        system.changeConstraint(1, changed);
    });


    CountingMaterial::destroyed = 0;
    {
        StructSystem owned(2);
        owned.addMaterial(std::make_unique<CountingMaterial>(1));
    }
    check_equal("destructor releases owned material",
                CountingMaterial::destroyed, 1);

    CountingMaterial::destroyed = 0;
    {
        StructSystem modern(4);
        auto material = std::make_unique<CountingMaterial>(3);
        auto* observer = material.get();
        modern.addMaterial(std::move(material));
        check_equal("unique ownership transferred", material == nullptr, 1);
        check_equal("owned material remains addressable",
                    modern.materials().at(3).get() == observer, 1);
    }
    check_equal("unique-owned material released", CountingMaterial::destroyed, 1);

    CountingMaterial::destroyed = 0;
    {
        StructSystem duplicates(6);
        duplicates.addMaterial(std::make_unique<CountingMaterial>(21));
        check_throws<std::invalid_argument>("reject duplicate material id", [&] {
            duplicates.addMaterial(std::make_unique<CountingMaterial>(21));
        });
        check_equal("duplicate material preserves original",
                    static_cast<int>(duplicates.materials().size()), 1);
        check_equal("rejected duplicate material released",
                    CountingMaterial::destroyed, 1);
    }
    check_equal("stored material released after scope", CountingMaterial::destroyed, 2);

    CountingMaterial::destroyed = 0;
    {
        StructSystem bulk(5);
        std::vector<std::unique_ptr<BaseMaterial>> materials;
        materials.push_back(std::make_unique<CountingMaterial>(10));
        materials.push_back(std::make_unique<CountingMaterial>(11));
        bulk.addMaterials(std::move(materials));
        check_equal("bulk unique materials stored",
                    static_cast<int>(bulk.materials().size()), 2);
    }
    check_equal("bulk unique materials released", CountingMaterial::destroyed, 2);

    const auto temp = std::filesystem::temp_directory_path() /
                      "project-nonlinear-modern-cpp-fs";
    std::filesystem::remove_all(temp);
    const auto nested = temp / "a" / "b";
    check_equal("create nested directory", makedir(nested.string()), 1);
    check_equal("detect directory", isDirExists(nested.string()), 1);

    const auto file = nested / "value.txt";
    {
        std::ofstream out(file);
        out << "vfife";
    }
    check_equal("detect regular file", isFileExists(file.string()), 1);
    check_equal("file is not directory", isDirExists(file.string()), 0);
    check_equal("directory is not file", isFileExists(nested.string()), 0);
    std::filesystem::remove_all(temp);

    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
