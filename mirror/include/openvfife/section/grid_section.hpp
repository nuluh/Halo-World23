// openVFIFE public API: experimental section-grid integration
#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "Eigen/Dense"
#include <openvfife/material/uniaxial.hpp>
#include <openvfife/numerics/quadrature.hpp>
#include <openvfife/section/section.hpp>

namespace openvfife {
// define a struct to restore point coordinate
struct SectionPoint { double y, z; };



struct CellNodes
{
    /* store cell corner nodes, WARNING: Please follow the instructions of node
    *  sequence. Again, FOLLOW THE SEQUENCE AS ILLUSTRAED IN THE FIGURE.
    *  ^
    *  |
    *  |    p1 _______ p2
    *  |       |     |
    *  |       | * * |
    *  |       |     |
    *  |       | * * |
    *  |    p3 |_____| p4
    *  |
    *  |------------> z
    */
    SectionPoint p1, p2, p3, p4;
};


struct SectionGridParameter
{
    /* Parameters for section grid generation
    * sections - number of sections
    * cells - number of cells
    * points - number of integration points of each cell(along one direction)
    * x1, x2, ..., x8 - paramers for grid generation for different section
    *
    * Retangle Section: only x1, x2 are used
    *   x1 - number of cells along section width
    *   x2 - number of cells along section height
    *
    * AngleSteel Section: only x1, x2 are used
    *   x1 - number of cells along vertical leg (including corner)
    *   x2 - number of cells along horizontal leg
    *
    * You can use x1, ..., x8 for other types of Section which is regular.
    * However, Please add a explanation here for other users.
    */
    int sections;       // number of sections
    int cells;          // number of cells
    int points;         // number of integration points of each cell
    int x1, x2, x3, x4, x5, x6, x7, x8; // control paramers
};


class SectionCell
{
    /*
    * only quadrilateral cell has been implemented here. If you want to use
    * other type of SectionCell such as triangle, please add new type of
    * SectionCell
    *
    *  ^ y
    *  |
    *  |    p1 _______ p2             * integration points
    *  |       |     |               -------
    *  |       | * * |               |     |  section cell
    *  |       |     |               |     |
    *  |       | * * |               |_____|
    *  |    p3 |_____| p4
    *  |
    *  |------------> z
    *
    * TODO:
    * Notes: The SectionCell class only support same material, however it is
    *        quite plain to assign a different type of material to each
    *        integration point. Just modified material_ and add a function for
    *        assign materials. and change SectionRespose::setCellMaterial()
    */
    private:
        int number_integral_points_{};// number of integral points along one axis
        int fractured_fiber{};
        const BaseMaterial* material_{nullptr};
        std::unique_ptr<GaussLegendre> legendre_;
        inline void setIntegratePointCoordinate();

    public:
        // record coordinates of integral points
        std::vector<std::vector<SectionPoint>> integral_points;
        SectionPoint cp1, cp2, cp3, cp4;   // cell node coordinates
        double area{};          // cell area

        // constitutive state at each integration point
        std::vector<std::vector<UniaxialMaterialState>> material_states;

        // restore degree of damage of cell
        std::vector<std::vector<bool>> damage;     // if fiber fractured or not
        double damage_degree{};

        // record total internal force, [N, miy, miz, mjy, mjz]
        std::vector<double> force;

    public:
        // points - integral points
        SectionCell(const BaseMaterial& mat, const CellNodes &cn, int p=2);
        SectionCell(const CellNodes &cn, int p=2);
        virtual ~SectionCell() = default;

        void setMaterial(const BaseMaterial& mat);
        virtual void response(const Eigen::RowVectorXd &transfer,
            const Eigen::VectorXd &deform);
};


class SectionResponse
{
    /* Section is composed by cells, see RectangleGrid. SectionResponse class
    * is responsible for gathering cells results.
    *
    * damageDegree() is trying to calculate the damage degree of the section,
    * which is defined as follow:
    *   $ \lambda = \frac{N_{collapse}}{N} $
    * $ N_{collapse} $ - collapsed integration points
    * $ N $ - all integration points
    *
    */
    public:
        std::vector<std::unique_ptr<SectionCell>> cells;
        std::vector<double> force;

    public:
        virtual ~SectionResponse() = default;
        virtual void setCellMaterial(const BaseMaterial& mat);
        virtual double damageDegree() const;
        virtual void response(const Eigen::RowVectorXd &transfer,
            const Eigen::VectorXd &deform);
};


class SectionGrid
{
    // interface of GridSection
    public:
        virtual ~SectionGrid() = default;
        virtual void initiateResponseContainers(const BaseSection& sect,
            SectionResponse& response) const = 0;
};


class RectangleGrid: public SectionGrid
{
    // only support same material type
    /* Grid Parameters
    * x1 - number of cells along vertical direction (y)
    * x2 - number of cells along horizontal direction (z)
    *
    *              ^ y
    *        ______|______
    *        |     |     |
    *        |  1  |  2  |
    *        |     |     |  height
    *    ----|-----------|-------> z
    *        |     |     |
    *        |  3  |  4  |
    *        |_____|_____|
    *              |
    *           width
    *
    */

    private:
        // rows - number of section cells along width
        // cols - number of section cells along height
        int rows_, cols_;
        int integrals_points_;

        inline void checkSectionType(const BaseSection& sect) const;

    public:
        RectangleGrid(const SectionGridParameter &p);
        ~RectangleGrid() override = default;

        void initiateResponseContainers(const BaseSection& sect,
            SectionResponse& response) const override;
};


class AngleSteelGrid: public SectionGrid
{
    /* Grid Parameters
    * x1 - number of cells along vertical direction (y), exclude the corner
    * x2 - number of cells along horizontal direction (z), exclude the corner
    *
    *         t1    ^ y
    *       _____   |
    *       |   |   |
    *       | 1 |   |    1______2
    *       |___|   |     |    |
    *  w1   |   |   |     |    |
    *       | 2 |   |    3|____|4
    *       |___|   |
    *       |   |   |
    * ------|-3-|---|---------------------> z
    *       |___|___|______________
    *       |   |   1  |  2  |  3  |
    *       |___|___|__|_____|_____|  t2
    *               |
    *               | w2
    *
    */
    private:
        int x1_, x2_;
        int integrals_points_;

        inline void checkSectionType(const BaseSection& sect) const;

    public:
        AngleSteelGrid(const SectionGridParameter &para);
        ~AngleSteelGrid() override = default;

        void initiateResponseContainers(const BaseSection& sect,
            SectionResponse& response) const override;
};


// creator of GridSection
class CreateGrid
{
    /* C++ has no reflection mechanism, which means the program can not create
    * the GridSection by the section type. CreateGridSection gives a very
    * simple implementation for it. If you create a new GridSetion class,
    * Please add your
    */
    public:
        static std::unique_ptr<SectionGrid> getInstance(const BaseSection& sect,
            const SectionGridParameter &para)
        {
            if (sect.type() == "Rectangle")
            {
                return std::make_unique<RectangleGrid>(para);
            }
            if (sect.type() == "AngleSteel")
            {
                return std::make_unique<AngleSteelGrid>(para);
            }
            throw std::invalid_argument("unsupported section type: " + sect.type());
        }
};

} // namespace openvfife
