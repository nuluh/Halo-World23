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

// $Date: 2020-05-22 $
// Written: Tan Biao
// Revised:
//
// Purpose: This file contains the class definition for Particle

// The interface:
//

#include <openvfife/numerics/quadrature.hpp>

#include <stdexcept>
#include <cmath>

using namespace std;

namespace openvfife {


/*****************************************************************************/
/*                             BaseIntegrator                                */
/*****************************************************************************/
BaseIntegrator::BaseIntegrator(int nodes)
    : nodes_(nodes)
{
    if (nodes <= 0)
    {
        throw std::invalid_argument("integration point count must be positive");
    }
}

double BaseIntegrator::x(int i) const
{
    if (i < 0 || static_cast<size_t>(i) >= x_.size())
    {
        throw std::out_of_range("integration point index out of range");
    }
    return x_[i];
}

double BaseIntegrator::w(int i) const
{
    if (i < 0 || static_cast<size_t>(i) >= w_.size())
    {
        throw std::out_of_range("integration weight index out of range");
    }
    return w_[i];
}

double BaseIntegrator::integrate1d(const vector<double> &val, double x1,
    double x2) const
{
    if (val.size() < static_cast<size_t>(nodes_))
    {
        throw std::invalid_argument("insufficient 1D quadrature values");
    }
    double sum = 0;
    for (int i = 0; i < nodes_; i++)
    {
        sum += val[i] * w_[i];
    }
    return (x2 - x1) / 2.0 * sum;
}

double BaseIntegrator::integrate2d(const vector<vector<double>> &val,
    double x1, double x2, double y1, double y2) const
{
    if (val.empty())
    {
        throw std::invalid_argument("empty 2D quadrature values");
    }
    if ((val.size() < static_cast<size_t>(nodes_)) ||
        (val[0].size() < static_cast<size_t>(nodes_)))
    {
        throw std::invalid_argument("insufficient 2D quadrature values");
    }
    double sum = 0;
    for (int i = 0; i < nodes_; i++)
    {
        for (int j = 0; j < nodes_; j++)
        {
            sum += w_[i] * w_[j] * val[i][j];
        }
    }
    return (x2 - x1) * (y2 - y1) / 4.0 * sum;
}

double BaseIntegrator::integrate3d(const vector<vector<vector<double>>> &val,
    double x1, double x2, double y1, double y2, double z1, double z2) const
{
    if (val.empty() || val[0].empty())
    {
        throw std::invalid_argument("empty 3D quadrature values");
    }
    if ((val.size() < static_cast<size_t>(nodes_)) ||
        (val[0].size() < static_cast<size_t>(nodes_)) ||
        (val[0][0].size() < static_cast<size_t>(nodes_)))
    {
        throw std::invalid_argument("insufficient 3D quadrature values");
    }
    double sum = 0;
    for (int i = 0; i < nodes_; i++)
    {
        for (int j = 0; j < nodes_; j++)
        {
            for (int k = 0; k < nodes_; k++)
            {
                sum += w_[i] * w_[j] * w_[k] * val[i][j][k];
            }
        }
    }
    return (x2 - x1) * (y2 - y1) * (z2 - z1) / 8.0 * sum;
}


/*****************************************************************************/
/*                              GaussLegendre                                */
/*****************************************************************************/
GaussLegendre::GaussLegendre(int nodes): BaseIntegrator(nodes)
{
    if (nodes_ > 10)
    {
        throw std::invalid_argument("Gauss-Legendre supports at most 10 points");
    }
    switch (nodes_)
    {
    case 1:
        x_.push_back(0.000000000000000);
        w_.push_back(2.000000000000000);
        break;
    case 2:
        x_.push_back(-0.577350269189625);
        x_.push_back(0.577350269189625);
        w_.push_back(1.000000000000000);
        w_.push_back(1.000000000000000);
        break;
    case 3:
        x_.push_back(-0.774596669241483);
        x_.push_back(0.000000000000000);
        x_.push_back(0.774596669241483);
        w_.push_back(0.555555555555555);
        w_.push_back(0.888888888888888);
        w_.push_back(0.555555555555555);
        break;
    case 4:
        x_.push_back(-0.861136311594052);
        x_.push_back(-0.339981043584856);
        x_.push_back(0.339981043584856);
        x_.push_back(0.861136311594052);
        w_.push_back(0.347854845137453);
        w_.push_back(0.652145154862546);
        w_.push_back(0.652145154862546);
        w_.push_back(0.347854845137453);
        break;
    case 5:
        x_.push_back(-0.906179845938663);
        x_.push_back(-0.538469310105683);
        x_.push_back(0.000000000000000);
        x_.push_back(0.538469310105683);
        x_.push_back(0.906179845938663);
        w_.push_back(0.236926885056189);
        w_.push_back(0.478628670499366);
        w_.push_back(0.568888888888888);
        w_.push_back(0.478628670499366);
        w_.push_back(0.236926885056189);
        break;
    case 6:
        x_.push_back(-0.932469514203152);
        x_.push_back(-0.661209386466264);
        x_.push_back(-0.238619186083196);
        x_.push_back(0.238619186083196);
        x_.push_back(0.661209386466264);
        x_.push_back(0.932469514203152);
        w_.push_back(0.171324492379170);
        w_.push_back(0.360761573048138);
        w_.push_back(0.467913934572691);
        w_.push_back(0.467913934572691);
        w_.push_back(0.360761573048138);
        w_.push_back(0.171324492379170);
        break;
    case 7:
        x_.push_back(-0.949107912342758);
        x_.push_back(-0.741531185599394);
        x_.push_back(-0.405845151377397);
        x_.push_back(0.000000000000000);
        x_.push_back(0.405845151377397);
        x_.push_back(0.741531185599394);
        x_.push_back(0.949107912342758);
        w_.push_back(0.129484966168869);
        w_.push_back(0.279705391489276);
        w_.push_back(0.381830050505118);
        w_.push_back(0.417959183673469);
        w_.push_back(0.381830050505118);
        w_.push_back(0.279705391489276);
        w_.push_back(0.129484966168869);
        break;
    case 8:
        x_.push_back(-0.960289856497536);
        x_.push_back(-0.796666477413626);
        x_.push_back(-0.525532409916328);
        x_.push_back(-0.183434642495649);
        x_.push_back(0.183434642495649);
        x_.push_back(0.525532409916328);
        x_.push_back(0.796666477413626);
        x_.push_back(0.960289856497536);
        w_.push_back(0.101228536290376);
        w_.push_back(0.222381034453374);
        w_.push_back(0.313706645877887);
        w_.push_back(0.362683783378361);
        w_.push_back(0.362683783378361);
        w_.push_back(0.313706645877887);
        w_.push_back(0.222381034453374);
        w_.push_back(0.101228536290376);
        break;
    case 9:
        x_.push_back(-0.968160239507626);
        x_.push_back(-0.836031107326635);
        x_.push_back(-0.613371432700590);
        x_.push_back(-0.324253423403808);
        x_.push_back(0.000000000000000);
        x_.push_back(0.324253423403808);
        x_.push_back(0.613371432700590);
        x_.push_back(0.836031107326635);
        x_.push_back(0.968160239507626);
        w_.push_back(0.081274388361574);
        w_.push_back(0.180648160694857);
        w_.push_back(0.260610696402935);
        w_.push_back(0.312347077040002);
        w_.push_back(0.330239355001259);
        w_.push_back(0.312347077040002);
        w_.push_back(0.260610696402935);
        w_.push_back(0.180648160694857);
        w_.push_back(0.081274388361574);
        break;
    case 10:
        x_.push_back(-0.973906528517171);
        x_.push_back(-0.865063366688984);
        x_.push_back(-0.679409568299024);
        x_.push_back(-0.433395394129247);
        x_.push_back(-0.148874338981631);
        x_.push_back(0.148874338981631);
        x_.push_back(0.433395394129247);
        x_.push_back(0.679409568299024);
        x_.push_back(0.865063366688984);
        x_.push_back(0.973906528517171);
        w_.push_back(0.066671344308688);
        w_.push_back(0.149451349150580);
        w_.push_back(0.219086362515982);
        w_.push_back(0.269266719309996);
        w_.push_back(0.295524224714752);
        w_.push_back(0.295524224714752);
        w_.push_back(0.269266719309996);
        w_.push_back(0.219086362515982);
        w_.push_back(0.149451349150580);
        w_.push_back(0.066671344308688);
        break;
    default:
        throw std::logic_error("unsupported Gauss-Legendre order");
    }
}


/*****************************************************************************/
/*                               GaussLobatto                                */
/*****************************************************************************/
GaussLobatto::GaussLobatto(int nodes): BaseIntegrator(nodes)
{
    if (nodes_ < 2 || nodes_ > 10)
    {
        throw std::invalid_argument("Gauss-Lobatto supports 2 to 10 points");
    }
    switch (nodes_)
    {
    case 2:
        x_.push_back(-1.000000000000000);
        x_.push_back(1.000000000000000);
        w_.push_back(1.000000000000000);
        w_.push_back(1.000000000000000);
        break;
    case 3:
        x_.push_back(-1.000000000000000);
        x_.push_back(0.000000000000000);
        x_.push_back(1.000000000000000);
        w_.push_back(0.333333333333333);
        w_.push_back(1.333333333333330);
        w_.push_back(0.333333333333333);
        break;
    case 4:
        x_.push_back(-1.000000000000000);
        x_.push_back(-0.447213595499957);
        x_.push_back(0.447213595499957);
        x_.push_back(1.000000000000000);
        w_.push_back(0.166666666666666);
        w_.push_back(0.833333333333333);
        w_.push_back(0.833333333333333);
        w_.push_back(0.166666666666666);
        break;
    case 5:
        x_.push_back(-1.000000000000000);
        x_.push_back(-0.654653670707977);
        x_.push_back(0.000000000000000);
        x_.push_back(0.654653670707977);
        x_.push_back(1.000000000000000);
        w_.push_back(0.100000000000000);
        w_.push_back(0.544444444444444);
        w_.push_back(0.711111111111111);
        w_.push_back(0.544444444444444);
        w_.push_back(0.100000000000000);
        break;
    case 6:
        x_.push_back(-1.000000000000000);
        x_.push_back(-0.765055323929464);
        x_.push_back(-0.285231516480645);
        x_.push_back(0.285231516480645);
        x_.push_back(0.765055323929464);
        x_.push_back(1.000000000000000);
        w_.push_back(0.066666666666667);
        w_.push_back(0.378474956297846);
        w_.push_back(0.554858377035486);
        w_.push_back(0.554858377035486);
        w_.push_back(0.378474956297846);
        w_.push_back(0.066666666666667);
        break;
    case 7:
        x_.push_back(-1.000000000000000);
        x_.push_back(-0.830223896278566);
        x_.push_back(-0.468848793470714);
        x_.push_back(0.000000000000000);
        x_.push_back(0.468848793470714);
        x_.push_back(0.830223896278566);
        x_.push_back(1.000000000000000);
        w_.push_back(0.047619047619048);
        w_.push_back(0.276826047361565);
        w_.push_back(0.431745381209862);
        w_.push_back(0.487619047619047);
        w_.push_back(0.431745381209862);
        w_.push_back(0.276826047361565);
        w_.push_back(0.047619047619048);
        break;
    case 8:
        x_.push_back(-1.000000000000000);
        x_.push_back(-0.871740148509606);
        x_.push_back(-0.591700181433142);
        x_.push_back(-0.209299217902478);
        x_.push_back(0.209299217902478);
        x_.push_back(0.591700181433142);
        x_.push_back(0.871740148509606);
        x_.push_back(1.000000000000000);
        w_.push_back(0.035714285714286);
        w_.push_back(0.210704227143506);
        w_.push_back(0.341122692483504);
        w_.push_back(0.412458794658703);
        w_.push_back(0.412458794658703);
        w_.push_back(0.341122692483504);
        w_.push_back(0.210704227143506);
        w_.push_back(0.035714285714286);
        break;
    case 9:
        x_.push_back(-1.000000000000000);
        x_.push_back(-0.899757995411460);
        x_.push_back(-0.677186279510737);
        x_.push_back(-0.363117463826178);
        x_.push_back(0.000000000000000);
        x_.push_back(0.363117463826178);
        x_.push_back(0.677186279510737);
        x_.push_back(0.899757995411460);
        x_.push_back(1.000000000000000);
        w_.push_back(0.027777777777778);
        w_.push_back(0.165495361560805);
        w_.push_back(0.274538712500161);
        w_.push_back(0.346428510973046);
        w_.push_back(0.371519274376417);
        w_.push_back(0.346428510973046);
        w_.push_back(0.274538712500161);
        w_.push_back(0.165495361560805);
        w_.push_back(0.027777777777778);
        break;
    case 10:
        x_.push_back(-1.000000000000000);
        x_.push_back(-0.919533908166458);
        x_.push_back(-0.738773865105505);
        x_.push_back(-0.477924949810444);
        x_.push_back(-0.165278957666387);
        x_.push_back(0.165278957666387);
        x_.push_back(0.477924949810444);
        x_.push_back(0.738773865105505);
        x_.push_back(0.919533908166458);
        x_.push_back(1.000000000000000);
        w_.push_back(0.022222222222222);
        w_.push_back(0.133305990851070);
        w_.push_back(0.224889342063126);
        w_.push_back(0.292042683679683);
        w_.push_back(0.327539761183897);
        w_.push_back(0.327539761183897);
        w_.push_back(0.292042683679683);
        w_.push_back(0.224889342063126);
        w_.push_back(0.133305990851070);
        w_.push_back(0.022222222222222);
        break;
    default:
        throw std::logic_error("unsupported Gauss-Lobatto order");
    }
}

} // namespace openvfife
