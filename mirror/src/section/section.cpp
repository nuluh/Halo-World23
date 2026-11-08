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

// $Date: 2020-05-11 $
// Written: Tan Biao
// Revised:
//
// Purpose: This file contains the class definition for Sections

// The interface:
//

#include <openvfife/section/section.hpp>
#include <openvfife/section/experimental.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <iostream>
using namespace std;

namespace openvfife {


/*****************************************************************************/
/*                                BaseSection                                */
/*****************************************************************************/
BaseSection::BaseSection(int id)
{
    type_ = "BaseSection";
    id_ = id;
    A_ = 0;
    Sy_ = 0;
    Sz_ = 0;
    CGy_ = 0;
    CGz_ = 0;
    SHy_ = 0;
    SHz_ = 0;
    Iyy_ = 0;
    Izz_ = 0;
    Iyz_ = 0;
    J_ = 0;
}

void BaseSection::initiate()
{
    calcArea();
    calcSy();
    calcSz();
    calcSHy();
    calcSHz();
    calcCGy();
    calcCGz();
    calcIyy();
    calcIzz();
    calcIyz();
}

void BaseSection::info() const
{
    cout << "Area: " << A_ << endl;
    cout << "Sy: " << Sy_ << ", Sz: " << Sz_ << endl;
    cout << "CGy: " << CGy_ << ", CGz: " << CGz_ << endl;
    cout << "SHy: " << SHy_ << ", SHz: " << SHz_ << endl;
    cout << "Iyy: " << Iyy_ << ", Izz: " << Izz_ << "Iyz: " << Iyz_ << endl;
}

int BaseSection::id() const
{
    return id_;
}

string BaseSection::type() const
{
    return type_;
}

double BaseSection::A() const
{
    return A_;
}

double BaseSection::Sy() const
{
    return Sy_;
}

double BaseSection::Sz() const
{
    return Sz_;
}

double BaseSection::CGy() const
{
    return CGy_;
}

double BaseSection::CGz() const
{
    return CGz_;
}

double BaseSection::SHy() const
{
    return SHy_;
}

double BaseSection::SHz() const
{
    return SHz_;
}

double BaseSection::Iyy() const
{
    return Iyy_;
}

double BaseSection::Izz() const
{
    return Izz_;
}

double BaseSection::Iyz() const
{
    return Iyz_;
}

double BaseSection::J() const
{
    return J_;
}


/*****************************************************************************/
/*                               CustomSection                               */
/*****************************************************************************/
CustomSection::CustomSection(int id, const CustomSectionParas &paras):
    BaseSection(id)
{
    type_ = "CustomSection";
    if (paras.A < 0.0)
    {
        throw std::invalid_argument("custom section area must be non-negative");
    }
    id_ = id;
    A_ = paras.A;
    Sy_ = paras.Sy;
    Sz_ = paras.Sz;
    SHy_ = paras.SHy;
    SHz_ = paras.SHz;
    CGy_ = paras.CGy;
    CGz_ = paras.CGz;
    Iyy_ = paras.Iyy;
    Izz_ = paras.Izz;
    Iyz_ = paras.Iyz;
    if (paras.J <= 0.0)
    {
        J_ = Iyy_ + Izz_;    // polar approximation when no J is supplied
    }
    else
    {
        J_ = paras.J;
    }
}

void CustomSection::calcArea()
{

}

void CustomSection::calcSy()
{

}

void CustomSection::calcSz()
{

}

void CustomSection::calcCGy()
{

}

void CustomSection::calcCGz()
{

}

void CustomSection::calcSHy()
{

}

void CustomSection::calcSHz()
{

}

void CustomSection::calcIyy()
{

}

void CustomSection::calcIzz()
{

}

void CustomSection::calcIyz()
{

}

void CustomSection::info() const
{
    cout << "Section ID: " << id_ << ", Type: " << type_ << endl;
    BaseSection::info();
}


/*****************************************************************************/
/*                                   Ctube                                   */
/*****************************************************************************/
Ctube::Ctube(int id, double ri, double ro):BaseSection(id)
{
    type_ = "Ctube";
    if ((ri > ro) || (ri < 0.0) || (ro <= 0.0))
    {
        throw invalid_argument("Please make sure ro > ri >= 0");
    }
    id_ = id;
    ri_ = ri;
    ro_ = ro;
    initiate();
}

void Ctube::calcArea()
{
    A_ = std::numbers::pi_v<double> * (ro_ * ro_ - ri_ * ri_);
}

void Ctube::calcSy()
{
    Sy_ = 0.0;
}

void Ctube::calcSz()
{
    Sz_ = 0.0;
}

void Ctube::calcSHy()
{
    cout << "calcSHy function have not defined!" << endl;
    SHy_ = 0.0;
}

void Ctube::calcSHz()
{
    cout << "calcSHz function have not defined!" << endl;
    SHz_ = 0.0;
}

void Ctube::calcCGy()
{
    CGy_ = Sy_ / A_;
}

void Ctube::calcCGz()
{
    CGz_ = Sz_ / A_;
}

void Ctube::calcIyy()
{
    Iyy_ = std::numbers::pi_v<double> * pow(2.0, 4.0) * (pow(ro_, 4.0) - pow(ri_, 4.0)) / 64.0;
}

void Ctube::calcIzz()
{
    Izz_ = std::numbers::pi_v<double> * pow(2.0, 4.0) * (pow(ro_, 4.0) - pow(ri_, 4.0)) / 64.0;
}

void Ctube::calcIyz()
{
    Iyz_ = 0.0;
}

double Ctube::ri() const
{
    return ri_;
}

double Ctube::ro() const
{
    return ro_;
}

void Ctube::info() const
{
    cout << "Section ID: " << id_ << ", Type: " << type_ << endl;
    cout << "Ri: " << ri_ << ", Ro: " << ro_ << endl;
    BaseSection::info();
}


/*****************************************************************************/
/*                                AngleSteel                                 */
/*****************************************************************************/
AngleSteel::AngleSteel(int id, double w1, double w2, double t1, double t2):
    BaseSection(id)
{
    type_ = "AngleSteel";
    if ((w1 <= 0.0) || (w2 <= 0.0) || (t1 <= 0.0) || (t2 <= 0.0))
    {
        throw std::invalid_argument("angle dimensions must be positive");
    }
    if (t1 > w2 || t2 > w1)
    {
        throw std::invalid_argument("angle thickness must not exceed its leg dimension");
    }
    id_ = id;
    w1_ = w1;
    w2_ = w2;
    t1_ = t1;
    t2_ = t2;
    initiate();
}

void AngleSteel::calcArea()
{
    A_ = (w1_ * t1_) + (w2_ * t2_) - (t1_ * t2_);
}

void AngleSteel::calcSy()
{
    double A1 = w1_ * w2_;
    double z1 = w2_ / 2.0;
    double A2 = -(w1_ - t1_) * (w2_ - t2_);
    double z2 = t1_ + (w2_ - t1_) / 2.0;
    Sy_ = z1 * A1 + z2 * A2;
}

void AngleSteel::calcSz()
{
    double A1 = w1_ * w2_;
    double y1 = w1_ / 2.0;
    double A2 = -(w1_ - t1_) * (w2_ - t2_);
    double y2 = t2_ + (w1_ - t2_) / 2.0;
    Sz_ = y1 * A1 + y2 * A2;

}

void AngleSteel::calcSHy()
{
    cout << "calcSHy function have not defined!" << endl;
    SHy_ = 0.0;
}

void AngleSteel::calcSHz()
{
    cout << "calcSHz function have not defined!" << endl;
    SHz_ = 0.0;
}

void AngleSteel::calcCGy()
{
    CGy_ = Sz_ / A_;
}

void AngleSteel::calcCGz()
{
    CGz_ = Sy_ / A_;
}

void AngleSteel::calcIyy()
{
    double zc1 = t1_ / 2.0;
    double A1 = w1_ * t1_;
    double Iyc1 = w1_ * pow(t1_, 3.0) / 12.0 + pow(CGz_-zc1, 2.0) * A1;

    double zc2 = (w2_ - t1_) / 2.0 + t1_;
    double A2 = (w2_ - t1_) * t2_;
    double Iyc2 = t2_ * pow(w2_-t1_, 3.0) / 12.0 + pow(CGz_-zc2, 2.0) * A2;
    Iyy_ = Iyc1 + Iyc2;
}

void AngleSteel::calcIzz()
{
    // ineritia moment to zc
    double yc1 = w1_ / 2.0;
    double A1 = w1_ * t1_;
    double Izc1 = t1_ * pow(w1_, 3.0) / 12.0 + pow(CGy_-yc1, 2.0) * A1;

    double yc2 = t2_ / 2.0;
    double A2 = (w2_ - t1_) * t2_;
    double Izc2 = (w2_ - t1_) * pow(t2_, 3.0) / 12.0 + pow(CGy_-yc2, 2.0) * A2;
    Izz_ = Izc1 + Izc2;
}

void AngleSteel::calcIyz()
{
    double yc1 = w1_ / 2.0;
    double zc1 = t1_ / 2.0;
    double A1 = w1_ * t1_;

    double yc2 = t2_ / 2.0;
    double zc2 = (w2_ - t1_) / 2.0 + t1_;
    double A2 = (w2_ - t1_) * t2_;

    double a1 = yc1 - CGy_;
    double b1 = zc1 - CGz_;
    double a2 = yc2 - CGy_;
    double b2 = zc2 - CGz_;
    double Iyz1 = a1 * b1 * A1;
    double Iyz2 = a2 * b2 * A2;
    Iyz_ = Iyz1 + Iyz2;
}

double AngleSteel::w1() const
{
    return w1_;
}

double AngleSteel::w2() const
{
    return w2_;
}

double AngleSteel::t1() const
{
    return t1_;
}

double AngleSteel::t2() const
{
    return t2_;
}

void AngleSteel::info() const
{
    cout << "Section ID: " << id_ << ", Type: " << type_ << endl;
    cout << "  w1: " << w1_ << ", w2: " << w2_
         << ", t1: " << t1_ << ", t2: " << t2_ << endl;
    BaseSection::info();
}


/*****************************************************************************/
/*                                 Rectangle                                 */
/*****************************************************************************/
Rectangle::Rectangle(int id, double width, double height):
    BaseSection(id)
{
    type_ = "Rectangle";
    if ((width <= 0.0) || (height <= 0.0))
    {
        throw invalid_argument("Please make sure width, height > 0");
    }

    width_ = width;
    height_ = height;
    initiate();
}

void Rectangle::calcArea()
{
    A_ = width_ * height_;
}

void Rectangle::calcSy()
{
    Sy_ = 0.0;
}

void Rectangle::calcSz()
{
    Sz_ = 0.0;
}

void Rectangle::calcSHy()
{
    cout << "calcSHy function have not defined!" << endl;
    SHy_ = 0.0;
}

void Rectangle::calcSHz()
{
    cout << "calcSHz function have not defined!" << endl;
    SHz_ = 0.0;
}

void Rectangle::calcCGy()
{
    CGy_ = Sy_ / A_;
}

void Rectangle::calcCGz()
{
    CGz_ = Sz_ / A_;
}

void Rectangle::calcIyy()
{
    Iyy_ = height_ * pow(width_, 3.0) / 12.0;
}

void Rectangle::calcIzz()
{
    Izz_ = width_ * pow(height_, 3.0) / 12.0;
}

void Rectangle::calcIyz()
{
    Iyz_ = 0.0;
}

double Rectangle::width() const
{
    return width_;
}

double Rectangle::height() const
{
    return height_;
}

void Rectangle::info() const
{
    cout << "Section ID: " << id_ << ", Type: " << type_ << endl;
    cout << "  width: " << width_ << ", height: " << height_ << endl;
    BaseSection::info();
}

} // namespace openvfife
