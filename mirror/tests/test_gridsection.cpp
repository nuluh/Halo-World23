#include <openvfife/section/grid_section.hpp>
#include <openvfife/section/experimental.hpp>

#include <cmath>
#include <cstdio>
#include <stdexcept>

using openvfife::LinearElastic;
using openvfife::CustomSectionParas;
using openvfife::CustomSection;
using openvfife::AngleSteel;
using openvfife::Rectangle;
using openvfife::CellNodes;
using openvfife::SectionGridParameter;
using openvfife::SectionCell;
using openvfife::SectionResponse;
using openvfife::RectangleGrid;
using openvfife::CreateGrid;

namespace {
int checks = 0;
int failures = 0;

class TensionCutoffMaterial final : public LinearElastic
{
public:
    TensionCutoffMaterial(int id, double E, double limit)
        : LinearElastic(id), limit_(limit)
    {
        setE(E);
    }

    bool isFractured(double strain) const override
    {
        return strain > limit_;
    }

private:
    double limit_;
};

void check_near(const char* label, double actual, double expected,
                double tolerance)
{
    ++checks;
    const bool ok = std::abs(actual - expected) <= tolerance;
    if (!ok) {
        ++failures;
    }
    std::printf("  %-34s got % .6e want % .6e %s\n",
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
    LinearElastic material(1);
    material.setE(200.0e9);

    const CellNodes nodes{
        .p1 = {.y = 0.5, .z = -0.5},
        .p2 = {.y = 0.5, .z = 0.5},
        .p3 = {.y = -0.5, .z = -0.5},
        .p4 = {.y = -0.5, .z = 0.5},
    };

    check_throws<std::invalid_argument>("reject zero cell integration points", [&] {
        SectionCell invalid(material, nodes, 0);
    });
    check_throws<std::invalid_argument>("reject negative cell integration points", [&] {
        SectionCell invalid(nodes, -1);
    });

    SectionCell cell(material, nodes, 2);

    Eigen::RowVectorXd transfer(5);
    transfer << 1.0, 0.0, 0.0, 0.0, 0.0;
    Eigen::VectorXd deformation(5);
    deformation << 1.0e-3, 0.0, 0.0, 0.0, 0.0;

    cell.response(transfer, deformation);

    check_near("axial resultant", cell.force[0], 200.0e6, 1.0e-6);
    check_near("moment y-i", cell.force[1], 0.0, 1.0e-9);
    check_near("moment z-i", cell.force[2], 0.0, 1.0e-9);
    check_near("moment y-j", cell.force[3], 0.0, 1.0e-9);
    check_near("moment z-j", cell.force[4], 0.0, 1.0e-9);
    check_near("damage degree", cell.damage_degree, 0.0, 0.0);
    check_near("cell state stress", cell.material_states[0][0].stress,
               200.0e6, 1.0e-6);
    check_near("cell state strain", cell.material_states[0][0].strain,
               1.0e-3, 1.0e-12);

    TensionCutoffMaterial cutoff(5, 200.0e9, 1.0e-4);
    SectionCell damaged_cell(cutoff, nodes, 2);
    Eigen::RowVectorXd bending_transfer(5);
    bending_transfer << 0.0, 0.0, 1.0, 0.0, 0.0;
    Eigen::VectorXd bending_deformation(5);
    bending_deformation << 0.0, 0.0, 1.0e-3, 0.0, 0.0;
    damaged_cell.response(bending_transfer, bending_deformation);
    check_near("half-fiber damage degree", damaged_cell.damage_degree, 0.5, 1.0e-12);

    Rectangle rectangle(2, 1.0, 1.0);
    SectionGridParameter grid_params{
        .sections = 1,
        .cells = 4,
        .points = 2,
        .x1 = 2,
        .x2 = 2,
        .x3 = 0,
        .x4 = 0,
        .x5 = 0,
        .x6 = 0,
        .x7 = 0,
        .x8 = 0,
    };
    auto grid = CreateGrid::getInstance(rectangle, grid_params);
    check_near("grid factory returned object", grid ? 1.0 : 0.0, 1.0, 0.0);
    SectionResponse response;
    grid->initiateResponseContainers(rectangle, response);
    response.setCellMaterial(material);
    response.response(transfer, deformation);

    check_near("grid cell count", static_cast<double>(response.cells.size()), 4.0, 0.0);
    check_near("grid axial resultant", response.force[0], 200.0e6, 1.0e-6);

    auto invalid_rectangle = grid_params;
    invalid_rectangle.x1 = 0;
    invalid_rectangle.cells = 0;
    check_throws<std::invalid_argument>("reject zero rectangle rows", [&] {
        RectangleGrid invalid(invalid_rectangle);
    });

    auto mismatched_rectangle = grid_params;
    mismatched_rectangle.cells = 3;
    check_throws<std::invalid_argument>("reject mismatched rectangle cells", [&] {
        RectangleGrid invalid(mismatched_rectangle);
    });

    AngleSteel angle(3, 1.0, 1.0, 0.1, 0.1);
    check_throws<std::invalid_argument>("reject wrong section for rectangle grid", [&] {
        SectionResponse invalid_response;
        grid->initiateResponseContainers(angle, invalid_response);
    });

    CustomSectionParas custom_params{};
    custom_params.A = 1.0;
    CustomSection custom(4, custom_params);
    check_throws<std::invalid_argument>("factory rejects unsupported section", [&] {
        CreateGrid::getInstance(custom, grid_params);
    });
    SectionCell deferred(nodes, 2);
    check_throws<std::logic_error>("deferred cell requires material", [&] {
        deferred.response(transfer, deformation);
    });
    deferred.setMaterial(material);
    deferred.response(transfer, deformation);
    check_near("deferred cell responds after assignment", deferred.force[0],
               200.0e6, 1.0e-6);

    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
