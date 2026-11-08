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
// Purpose: This file contains the class definition for materials, including
// base class BaseMaterial (abstract base class), UniaxialPlastic, NdPlastic,
// and several derived materials - LinearElastic and uniaxial plastic models.
// UniBilinearStatic, UniBilinear, UniPower.

// The interface:
//

#include <openvfife/material/uniaxial.hpp>
#include <openvfife/material/experimental.hpp>

#include <stdexcept>
#include <cmath>
#include <iostream>
using namespace std;

namespace openvfife {

/*****************************************************************************/
/*                                  Helper                                   */
/*****************************************************************************/
inline int sign(int x)
{
    return (x > 0) - (x < 0);
}

inline int sign(double x)
{
    return (x > 0) - (x < 0);
}

/*****************************************************************************/
/*                               BaseMaterial                                */
/*****************************************************************************/
BaseMaterial::BaseMaterial(int id)
{
    id_ = id;
    type_ = "BaseMaterial";
    E_ = 0.0;
    G_ = 0.0;
    nu_ = 0.0;
    density_ = 0.0;
    ce_ = 0.0;
    damp_ratio_ = 0.0;
    damp_coef_ = 0.0;
    mu_ = 0.0;
    limited_tensile_strain_ = 9.99e9;
}

int BaseMaterial::id() const
{
    return id_;
}

string BaseMaterial::type() const
{
    return type_;
}

double BaseMaterial::E() const
{
    return E_;
}

double BaseMaterial::G() const
{
    return G_;
}

double BaseMaterial::nu() const
{
    return nu_;
}

double BaseMaterial::density() const
{
    return density_;
}

double BaseMaterial::ce() const
{
    return ce_;
}

double BaseMaterial::damp_ratio() const
{
    return damp_ratio_;
}

double BaseMaterial::damp_coef() const
{
    return damp_coef_;
}

double BaseMaterial::mu() const
{
    return mu_;
}

double BaseMaterial::yield_stress() const
{
    return yield_stress_;
}

double BaseMaterial::limited_tensile_strain() const
{
    return limited_tensile_strain_;
}

bool BaseMaterial::isFractured(double strain) const
{
    bool flag = false;
    if (strain > limited_tensile_strain_)
    {
        flag = true;
    }
    return flag;
}

void BaseMaterial::setE(double val)
{
    if (val < 0)
    {
        throw std::invalid_argument("Young's modulus must be non-negative");
    }
    E_ = val;
}

void BaseMaterial::setG(double val)
{
    if (val < 0)
    {
        throw std::invalid_argument("shear modulus must be non-negative");
    }
    G_ = val;
}

void BaseMaterial::setNu(double val)
{
    if (val < 0)
    {
        throw std::invalid_argument("Poisson ratio must be non-negative");
    }
    nu_ = val;
    G_ = E_ / 2.0 / (1 + nu_);
}

void BaseMaterial::setDensity(double val)
{
    if (val < 0)
    {
        throw std::invalid_argument("density must be non-negative");
    }
    density_ = val;
}

double BaseMaterial::calcCe() const
{
    // ce_ = sqrt(E_ / density_);
    // return ce_;
    return sqrt(E_ / density_);
}

void BaseMaterial::setDampRatio(double val)
{
    damp_ratio_ = val;
}

void BaseMaterial::setDampCoef(double val)
{
    damp_coef_ = val;
}

void BaseMaterial::setMu(double val)
{
    mu_ = val;
}

void BaseMaterial::setLimitedTensileStrain(double tensile)
{
    if (tensile <= 0.0)
    {
        throw std::invalid_argument("limited tensile strain must be positive");
    }
    limited_tensile_strain_ = tensile;
}

void BaseMaterial::info() const
{
    cout << "Material ID: " << id_ << endl;
    cout << "Material Type: " << type_ << endl;
    cout << "E: " << E_ << endl;
    cout << "G: " << G_ << endl;
    cout << "density: " << density_ << endl;
    cout << "Ultimate tensile strain: " << limited_tensile_strain_ << endl;
}


/*****************************************************************************/
/*                              LinearElastic                                */
/*****************************************************************************/
LinearElastic::LinearElastic(int id): BaseMaterial(id)
{
    type_ = "LinearElastic";
}

double LinearElastic::update(double delta_strain,
    UniaxialMaterialState &state) const
{
    state.strain += delta_strain;
    state.stress += E_ * delta_strain;
    return E_;
}


/*****************************************************************************/
/*                              UniaxialPlastic                              */
/*****************************************************************************/
UniaxialPlastic::UniaxialPlastic(int id): BaseMaterial(id)
{
    type_ = "UniaxialPlastic";
    m_ = 0.0;
    yield_stress_ = 0.0;

}

void UniaxialPlastic::setYieldStress(double val)
{
    if (val <= 0)
    {
        throw std::invalid_argument("yield stress must be positive");
    }
    yield_stress_ = val;
}

void UniaxialPlastic::info() const
{
    BaseMaterial::info();
    cout << "yield stress: " << yield_stress_ << endl;
}


/*****************************************************************************/
/*                                 UniIdeal                                  */
/*****************************************************************************/
UniIdeal::UniIdeal(int id, double E, double yield_stress):UniaxialPlastic(id)
{
    if (E <= 0)
    {
        throw std::invalid_argument("UniIdeal requires positive Young's modulus");
    }
    if (yield_stress <= 0)
    {
        throw std::invalid_argument("UniIdeal requires positive yield stress");
    }
    type_ = "UniIdeal";
    E_ = E;
    Et_ = 0;
    yield_stress_ = yield_stress;
}

double UniIdeal::update(double delta_strain,
    UniaxialMaterialState &state) const
{
    /*  Return Mapping Algorithm, << Computational Inelasticity>> Simo */
    // step 1: database at last step, ${{\epsilo}_{n}^{p}, {\alpha}_n, q_n}$
    double Eq = E_;

    // step 2: strain field at current step
    state.strain += delta_strain;

    // step 3: compute elastic trial stress and test for plastic loading
    double stress_trial = E_ * (state.strain - state.plastic_strain);
    double yield_func_trial = fabs(stress_trial) - yield_stress_;

    // step 4: return-mapping
    if (yield_func_trial < 0)
    {
        // not yield
        state.stress = stress_trial;
    }
    else
    {
        // yield
        double dgamma = yield_func_trial / E_;
        state.stress = stress_trial - dgamma * E_ * sign(stress_trial);
        state.plastic_strain += dgamma * sign(stress_trial);
        Eq = Et_;
    }
    return Eq;
}


/*****************************************************************************/
/*                               UniBilinear                                 */
/*****************************************************************************/
UniBilinear::UniBilinear(int id, double E, double Et, double yield_stress,
    double m):UniaxialPlastic(id)
{
    if (E <= 0)
    {
        throw std::invalid_argument("UniBilinear requires positive Young's modulus");
    }
    if (Et < 0 || Et > E)
    {
        throw std::invalid_argument("UniBilinear tangent modulus must be in [0, E]");
    }
    if (yield_stress <= 0)
    {
        throw std::invalid_argument("UniBilinear requires positive yield stress");
    }
    if (m < 0 || m > 1)
    {
        throw std::invalid_argument("UniBilinear hardening mix must be in [0, 1]");
    }
    type_ = "UniBilinear";
    E_ = E;
    Et_ = Et;
    yield_stress_ = yield_stress;
    m_ = m;
}

double UniBilinear::Et() const
{
    return Et_;
}

void UniBilinear::info() const
{
    UniaxialPlastic::info();
    cout << "Et: " << Et_ << endl;
}

double UniBilinear::update(double delta_strain,
    UniaxialMaterialState &state) const
{
    /*  Return Mapping Algorithm, << Computational Inelasticity>> Simo
    */
    // step 1: database at last step, ${{\epsilo}_{n}^{p}, {\alpha}_n, q_n}$
    double Eq = E_;

    // step 2: strain field at current step
    state.strain += delta_strain;

    // step 3: compute elastic trial stress and test for plastic loading
    double stress_trial = E_ * (state.strain - state.plastic_strain);
    double kesi_trial = stress_trial - (1 - m_) * state.back_stress;
    double yield_func_trial = fabs(kesi_trial) - (yield_stress_ + m_ * state.hardening);

    // step 4: return-mapping
    if (yield_func_trial < 0)
    {
        // not yield
        state.stress = stress_trial;
    }
    else
    {
        // yield
        double H = E_ * Et_ / (E_ - Et_);
        double dgamma = yield_func_trial / (E_ + m_ * H + (1 - m_) * H);
        state.stress = stress_trial - dgamma * E_ * sign(kesi_trial);
        state.plastic_strain += dgamma * sign(kesi_trial);
        state.back_stress += dgamma * H * sign(kesi_trial);
        state.hardening += dgamma * H;
        Eq = Et_;
    }
    return Eq;
}


/*****************************************************************************/
/*                                 UniPower                                  */
/*****************************************************************************/
UniPower::UniPower(int id, double E, double n, double yield_stress, double m):
    UniaxialPlastic(id)
{
    if (E <= 0)
    {
        throw std::invalid_argument("UniPower requires positive Young's modulus");
    }
    if (yield_stress <= 0)
    {
        throw std::invalid_argument("UniPower requires positive yield stress");
    }
    if (m < 0 || m > 1)
    {
        throw std::invalid_argument("UniPower hardening mix must be in [0, 1]");
    }
    type_ = "UniPower";
    E_ = E;
    n_ = n;
    k_ = yield_stress / pow(yield_stress / E, n);
    yield_stress_ = yield_stress;
    m_ = m;
}

void UniPower::info() const
{
    UniaxialPlastic::info();
    cout << "k: " << k_ << ", n: " << n_ << endl;
}

double UniPower::update(double delta_strain,
    UniaxialMaterialState &state) const
{
    /*  Return Mapping Algorithm, << Computational Inelasticity>> Simo
    */
    // step 1: database at last step, ${{\epsilo}_{n}^{p}, {\alpha}_n, q_n}$
    double Eq = E_;

    // step 2: strain field at current step
    state.strain += delta_strain;

    // step 3: compute elastic trial stress and test for plastic loading
    double stress_trial = E_ * (state.strain - state.plastic_strain);
    double kesi_trial = stress_trial - (1 - m_) * state.back_stress;
    double yield_func_trial = fabs(kesi_trial) - (yield_stress_ + m_ * state.hardening);


    // step 4: return-mapping
    if (yield_func_trial < 0)
    {
        // not yield
        state.stress = stress_trial;
    }
    else
    {
        // yield
        double yield_strain = yield_stress_ / E_;
        // if (fabs(kesi_trial) < 0)
        // {
        //     double sigma = (1 - m_) * (*alpha) - (yield_stress_ + m_ * (*q));
        //     yield_strain = pow(sigma / k_, 1 / n_);
        // }
        // else
        // {
        //     double sigma = (1 - m_) * (*alpha) + (yield_stress_ + m_ * (*q));
        //     yield_strain = pow(sigma / k_, 1 / n_);
        // }
        double Ep = k_ * n_ * pow(fabs(state.strain - yield_strain), n_-1);
        double Et = E_ * Ep / (E_ + Ep);
        double dgamma = yield_func_trial / (E_ + m_ * Ep + (1 - m_) * Ep);
        state.stress += Et * delta_strain;
        // cout << "Ep * dgamma = " << Ep * dgamma * sign(kesi_trial);
        // cout << ", Et * ds = " << Et * ds << ", E * es = ";
        // cout << E_ * (ds - dgamma * sign(kesi_trial)) << endl;
        state.plastic_strain += dgamma * sign(kesi_trial);
        state.back_stress += dgamma * Ep * sign(kesi_trial);
        state.hardening += dgamma * Ep;
        Eq = Et;
        // Eq = H;
    }
    return Eq;
}


/*****************************************************************************/
/*                                 NdPlastic                                 */
/*****************************************************************************/
// NdPlastic::NdPlastic(int id): BaseMaterial(id)
// {
//     type_ = "NdPlastic";
//     Et_.setZero();
//     Eq_.setZero();
//     stress_.setZero();
//     strain_.setZero();
//     delta_strain_.setZero();
// }

// NdPlastic::~NdPlastic()
// {

// }

// void NdPlastic::info() const
// {
//     cout << "Material ID: " << id_ << endl;
//     cout << "Material Type: " << type_ << endl;
//     cout << "E: " << E_ << endl;
//     cout << "G: " << G_ << endl;
//     cout << "Et: " << Et_ <<endl;
//     cout << "density: " << density_ << endl;
// }

} // namespace openvfife
