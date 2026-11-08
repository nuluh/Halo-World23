/* ************************************************************************* **
**    OpenVFIFE - Open System for Vector Form Instrinsic                     **
**                Finite Element Method (VFIFE)                              **
**                GinkGo(Tan Biao)                                           **
**                                                                           **
**                                                                           **
** (C) Copyright 2021, The GinkGo(Tan Biao). All Rights Reserved.            **
**                                                                           **
** Commercial use of this program without express permission of              **
** GinkGo(Tan Biao), is strictly prohibited.  See                            **
** file 'COPYRIGHT'  in main directory for information on usage and          **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.                  **
**                                                                           **
** Developed by:                                                             **
**      Tan Biao (ginkgoltd@outlook.com)                                     **
**                                                                           **
** ************************************************************************* */

// $Date: 2021-06-07$
// Written: Tan Biao
// Revised:
//
// Purpose: This file contains the class definition for SectionGrid, including
// SectionCell, SectionResponse, RectangleGrid, AngleSteelGrid

// The interface:
//

#include <openvfife/section/grid_section.hpp>
#include <openvfife/section/experimental.hpp>

#include <stdexcept>


using namespace std;

namespace openvfife {


/*****************************************************************************/
/*                                  Helpers                                  */
/*****************************************************************************/
inline void initCellDamage(vector<vector<bool>> &v, int rows, int cols)
{
    v.assign(static_cast<size_t>(rows),
             vector<bool>(static_cast<size_t>(cols), false));
}


/*****************************************************************************/
/*                                 SectionCell                               */
/*****************************************************************************/
inline void SectionCell::setIntegratePointCoordinate()
{
    // initiate coordiantes of integral points
    // [[(y1, z1), (y1, z2),... (y1, zn)],
    //  [(y2, z1), (y2, z2),... (y2, zn)],
    //  ...
    //  [(yn, z1), (yn, z2),... (yn, zn)]]
    for (int i = 0; i < number_integral_points_; i++)
    {
        vector<SectionPoint> tmp;
        for (int j = 0; j < number_integral_points_; j++)
        {
            double y = (cp1.y - cp3.y) / 2.0 * legendre_->x(i) +
                       (cp1.y + cp3.y) / 2.0;
            double z = (cp2.z - cp1.z) / 2.0 * legendre_->x(j) +
                       (cp1.z + cp2.z) / 2.0;
            SectionPoint p = {.y = y, .z=z};
            tmp.push_back(p);
        }
        integral_points.push_back(tmp);
    }
}

SectionCell::SectionCell(const BaseMaterial& mat, const CellNodes &cn, int p)
{
    // set material
    material_ = &mat;

    // initiate integrator
    if (p <= 0)
    {
        throw std::invalid_argument("SectionCell integration points must be positive");
    }
    number_integral_points_ = p;
    legendre_ = std::make_unique<GaussLegendre>(p);

    // initiate cell nodes
    cp1 = cn.p1;
    cp2 = cn.p2;
    cp3 = cn.p3;
    cp4 = cn.p4;
    area = fabs((cn.p1.y - cn.p3.y) * (cn.p1.z - cn.p2.z));

    // initiate integral_points, after setting number_integral_points_
    setIntegratePointCoordinate();

    // initiate constitutive state at each integration point
    material_states.assign(static_cast<size_t>(p),
                           vector<UniaxialMaterialState>(static_cast<size_t>(p)));
    initCellDamage(damage, p, p);
    fractured_fiber = 0;

    // initiate force container
    for (int i = 0; i < 5; i++)
    {
        force.push_back(0);
    }
}

SectionCell::SectionCell(const CellNodes &cn, int p)
{
    // initiate integrator
    if (p <= 0)
    {
        throw std::invalid_argument("SectionCell integration points must be positive");
    }
    number_integral_points_ = p;
    legendre_ = std::make_unique<GaussLegendre>(p);

    // initiate cell nodes
    cp1 = cn.p1;
    cp2 = cn.p2;
    cp3 = cn.p3;
    cp4 = cn.p4;
    area = fabs((cn.p1.y - cn.p3.y) * (cn.p1.z - cn.p2.z));

    // initiate integral_points, after setting number_integral_points_
    setIntegratePointCoordinate();

    // initiate constitutive state at each integration point
    material_states.assign(static_cast<size_t>(p),
                           vector<UniaxialMaterialState>(static_cast<size_t>(p)));
    initCellDamage(damage, p, p);
    fractured_fiber = 0;

    // initiate force container
    for (int i = 0; i < 5; i++)
    {
        force.push_back(0);
    }
}

void SectionCell::setMaterial(const BaseMaterial& mat)
{
    material_ = &mat;
}

void SectionCell::response(const Eigen::RowVectorXd &transfer,
    const Eigen::VectorXd &deform)
{
    if (material_ == nullptr)
    {
        throw std::logic_error("SectionCell material must be assigned before response");
    }

    // step 1: calculate strain increment and response of each integral points
    for (int i = 0; i < number_integral_points_; i++)
    {
        for (int j = 0; j < number_integral_points_; j++)
        {
            // compute the strain at each integral points
            double ds = transfer(0) * deform(0)
                      + transfer(1) * deform(1) * integral_points[i][j].z
                      + transfer(2) * deform(2) * integral_points[i][j].y
                      + transfer(3) * deform(3) * integral_points[i][j].z
                      + transfer(4) * deform(4) * integral_points[i][j].y;

            auto& state = material_states[i][j];
            if (damage[i][j])
            {
                // fiber already fractured
                state.stress = 0;
            }
            else
            {
                // fiber still working
                if (material_->isFractured(state.strain + ds))
                {
                    // fiber fracture
                    damage[i][j] = true;
                    state.stress = 0;
                    fractured_fiber += 1;
                }
                else
                {
                    // get material response of each integration point
                    material_->update(ds, state);
                }
            }
        }
    }
    const double total_fibers = static_cast<double>(number_integral_points_)
                                * number_integral_points_;
    damage_degree = static_cast<double>(fractured_fiber) / total_fibers;

    // step 2: integration of cell stress
    // reset SectionResponse.force
    for (int i = 0; i < 5; i++)
    {
        force[i] = 0;   // fjx, miy, miz, mjy, mjz
    }

    double factor = (cp1.y - cp3.y) * (cp2.z - cp1.z) / 4.0;
    for (int i = 0; i < number_integral_points_; i++)
    {
        for (int j = 0; j < number_integral_points_; j++)
        {
            force[0] += legendre_->w(i) * legendre_->w(j) * transfer(0)
                        * factor * material_states[i][j].stress;
            force[1] += legendre_->w(i) * legendre_->w(j) * transfer(1)
                        * integral_points[i][j].z * factor * material_states[i][j].stress;
            force[2] += legendre_->w(i) * legendre_->w(j) * transfer(2)
                        * integral_points[i][j].y * factor * material_states[i][j].stress;
            force[3] += legendre_->w(i) * legendre_->w(j) * transfer(3)
                        * integral_points[i][j].z * factor * material_states[i][j].stress;
            force[4] += legendre_->w(i) * legendre_->w(j) * transfer(4)
                        * integral_points[i][j].y * factor * material_states[i][j].stress;
        }
    }
}


/*****************************************************************************/
/*                               SectionResponse                             */
/*****************************************************************************/
void SectionResponse::setCellMaterial(const BaseMaterial& mat)
{
    for (auto& cell : cells)
    {
        cell->setMaterial(mat);
    }
}

double SectionResponse::damageDegree() const
{
    double damage = 0;
    for (const auto& cell : cells)
    {
        damage += cell->damage_degree;
    }
    return damage / cells.size();
}

void SectionResponse::response(const Eigen::RowVectorXd &transfer,
    const Eigen::VectorXd &deform)
{
    // step 1: calculate cell response
    for (auto& cell : cells)
    {
        cell->response(transfer, deform);
    }

    for (auto& value : force)
    {
        value = 0;
    }

    for (std::size_t i = 0; i < force.size(); ++i)
    {
        for (const auto& cell : cells)
        {
            force[i] += cell->force[i];
        }
    }
}


/*****************************************************************************/
/*                                 GridSection                               */
/*****************************************************************************/
// interface


/*****************************************************************************/
/*                               GridRectangle                               */
/*****************************************************************************/
RectangleGrid::RectangleGrid(const SectionGridParameter &p): SectionGrid()
{
    if (p.x1 <= 0 || p.x2 <= 0)
    {
        throw std::invalid_argument("rectangle grid dimensions must be positive");
    }
    if (p.points <= 0)
    {
        throw std::invalid_argument("rectangle grid integration points must be positive");
    }
    if (p.x1 * p.x2 != p.cells)
    {
        throw std::invalid_argument("rectangle grid cell count does not match dimensions");
    }

    rows_ = p.x1;
    cols_ = p.x2;
    integrals_points_ = p.points;
}

inline void RectangleGrid::checkSectionType(const BaseSection& sect) const
{
    // check section type
    if (sect.type() != "Rectangle")
    {
        throw std::invalid_argument("RectangleGrid requires a Rectangle section");
    }
}

void RectangleGrid::initiateResponseContainers(const BaseSection& sect,
    SectionResponse& response) const
{
    // check section type first
    checkSectionType(sect);

    // coordinate of cell nodes
    const auto& rectangle = static_cast<const Rectangle&>(sect);
    const double w = rectangle.width();
    const double h = rectangle.height();

    vector<vector<SectionPoint>> sp;
    for (int i = 0; i <= rows_; i++)
    {
        vector<SectionPoint> tmp;
        for (int j = 0; j <= cols_; j++)
        {
            double y = h / 2.0 - h / cols_ * i;
            double z = -w / 2.0 + w / rows_ * j;
            SectionPoint p = {.y=y, .z=z};
            tmp.push_back(p);
        }
        sp.push_back(tmp);
    }

    // create cells
    for (int i = 0; i < rows_; i++)
    {
        for (int j = 0; j < cols_; j++)
        {
            CellNodes cn = {.p1=sp[i][j], .p2=sp[i][j+1], .p3=sp[i+1][j],
                            .p4=sp[i+1][j+1]};
            response.cells.push_back(std::make_unique<SectionCell>(cn, integrals_points_));
        }
    }

    // initiate force container
    for (int i = 0; i < 5; i++)
    {
        response.force.push_back(0);
    }
}


/*****************************************************************************/
/*                              GridAngleSteel                               */
/*****************************************************************************/
AngleSteelGrid::AngleSteelGrid(const SectionGridParameter &para):
    SectionGrid()
{
    if (para.x1 <= 0 || para.x2 <= 0)
    {
        throw std::invalid_argument("angle grid dimensions must be positive");
    }
    if (para.points <= 0)
    {
        throw std::invalid_argument("angle grid integration points must be positive");
    }
    if ((para.x1 + para.x2) != para.cells)
    {
        throw std::invalid_argument("angle grid cell count does not match dimensions");
    }

    x1_ = para.x1;
    x2_ = para.x2;
    integrals_points_ = para.points;
}

inline void AngleSteelGrid::checkSectionType(const BaseSection& sect) const
{
    // check section type
    if (sect.type() != "AngleSteel")
    {
        throw std::invalid_argument("AngleSteelGrid requires an AngleSteel section");
    }
}

void AngleSteelGrid::initiateResponseContainers(const BaseSection& sect,
    SectionResponse& response) const
{
    // check section type first
    checkSectionType(sect);

    // coordinate of cell nodes, multiple 10e6 for eliminate truncation error
    const auto& angle = static_cast<const AngleSteel&>(sect);
    const double w1 = angle.w1();
    const double w2 = angle.w2();
    const double t1 = angle.t1();
    const double t2 = angle.t2();
    const double yc = angle.CGy();
    const double zc = angle.CGz();

    // vertical part
    for (int i = 0; i < x1_; i++)
    {
        // point 1
        double y1 = w1 - yc - (w1 - t2) / x1_ * i;
        double z1 = -zc ;
        SectionPoint p1 = {.y=y1, .z=z1};

        // point 2
        double y2 = y1;
        double z2 = z1 + t1;
        SectionPoint p2 = {.y=y2, .z=z2};

        // point 3
        double y3 = y1 - (w1 - t2) / x1_;
        double z3 = z1;
        SectionPoint p3 = {.y=y3, .z=z3};

        // point 4
        double y4 = y3;
        double z4 = z2;
        SectionPoint p4 = {.y=y4, .z=z4};

        // create cell
        CellNodes cn = {.p1=p1, .p2=p2, .p3=p3, .p4=p4};
        response.cells.push_back(std::make_unique<SectionCell>(cn, integrals_points_));
    }

    // corner part
    {
        // point 1
        double y1 = t2 - yc;
        double z1 = -zc;
        SectionPoint p1 = {.y=y1, .z=z1};

        // point 2
        double y2 = y1;
        double z2 = z1 + t1;
        SectionPoint p2 = {.y=y2, .z=z2};

        // point 3
        double y3 = y1 - t2;
        double z3 = z1;
        SectionPoint p3 = {.y=y3, .z=z3};

        // point 4
        double y4 = y3;
        double z4 = z2;
        SectionPoint p4 = {.y=y4, .z=z4};

        // create cell
        CellNodes cn = {.p1=p1, .p2=p2, .p3=p3, .p4=p4};
        response.cells.push_back(std::make_unique<SectionCell>(cn, integrals_points_));
    }

    // horizontal part
    for (int i = 0; i < x2_; i++)
    {
        // point 1
        double y1 = t2 - yc;
        double z1 = t1 - zc + (w2 - t1) / x2_ * i;
        SectionPoint p1 = {.y=y1, .z=z1};

        // point 2
        double y2 = y1;
        double z2 = z1 + (w2 - t1) / x2_;
        SectionPoint p2 = {.y=y2, .z=z2};

        // point 3
        double y3 = y1 - t2;
        double z3 = z1;
        SectionPoint p3 = {.y=y3, .z=z3};

        // point 4
        double y4 = y3;
        double z4 = z2;
        SectionPoint p4 = {.y=y4, .z=z4};

        // create cell
        CellNodes cn = {.p1=p1, .p2=p2, .p3=p3, .p4=p4};
        response.cells.push_back(std::make_unique<SectionCell>(cn, integrals_points_));
    }

    // initiate force container
    for (int i = 0; i < 5; i++)
    {
        response.force.push_back(0);
    }
}

} // namespace openvfife
