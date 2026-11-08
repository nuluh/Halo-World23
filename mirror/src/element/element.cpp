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
// Purpose: This file contains the class definition for elements, including
// base class BaseElement (abstract base class), StructElement, PlaneElement,
// SolidElement, and several elements derived from StructElement - Link2D,
// Link2DLD, Link3D, Link3DLD, Beam2D and Beam3D.

// The interface:
//

#include <openvfife/element/base.hpp>
#include <openvfife/element/frame.hpp>
#include <openvfife/element/experimental.hpp>

#include <stdexcept>
#include <cmath>
#include <iostream>
using namespace std;

namespace openvfife {


/*****************************************************************************/
/*                                 helper                                    */
/*****************************************************************************/
inline double euclideanDistance(const StdArray6d &pi, const StdArray6d &pj)
{
    double ans = 0.0;
    for (int i = 0; i < 3; i++)
    {
        ans += (pi[i] - pj[i]) * (pi[i] - pj[i]);
    }
    return sqrt(ans);
};

inline void StdArray6dMulNum(StdArray6d &arr, double num)
{
    for (int i = 0; i < 6; i++)
    {
        arr[i] *= num;
    }
}

inline void StdAarrayAddNum(StdArray6d &arr, double num)
{
    for (int i = 0; i < 6; i++)
    {
        arr[i] += num;
    }
}


/*****************************************************************************/
/*                               BaseElement                                 */
/*****************************************************************************/
BaseElement::BaseElement(int id)
    : id_(id)
{
}

void BaseElement::addParticle(Particle* p)
{
    if (p == nullptr)
    {
        throw std::invalid_argument("element particle must not be null");
    }
    const auto required = static_cast<std::size_t>(num_of_particles_);
    if (particles_.size() >= required)
    {
        throw std::invalid_argument("element already has the required particle count");
    }
    if (particles_.contains(p->id()))
    {
        throw std::invalid_argument("element particle ids must be unique");
    }
    particles_.emplace(p->id(), p);
}

void BaseElement::info() const
{
    cout << "Element ID: " << id_ << endl;
    cout << "Element Type: " << type_ << endl;
    cout << "Particle ID: (";
    for (auto iter=particles_.begin(); iter != particles_.end(); iter++)
    {
        cout << iter->first << ", ";
    }
    cout << ")" << endl;
    cout << "\nMaterial: " << endl;
    material_->info();
    cout << "\nSection: " << endl;
    section_->info();
    cout << "\n" << endl;
}

int BaseElement::id() const
{
    return id_;
}

string BaseElement::type() const
{
    return type_;
}

int BaseElement::num_of_particles() const
{
    return num_of_particles_;
}

const map<int, Particle*>& BaseElement::particles() const noexcept
{
    return particles_;
}

double BaseElement::mass() const
{
    return mass_;
}

const BaseMaterial& BaseElement::material() const noexcept
{
    return *material_;
}

const BaseSection& BaseElement::section() const noexcept
{
    return *section_;
}

const Eigen::Vector3d& BaseElement::ex() const noexcept
{
    return ex_;
}

const Eigen::Vector3d& BaseElement::ey() const noexcept
{
    return ey_;
}

const Eigen::Vector3d& BaseElement::ez() const noexcept
{
    return ez_;
}

const Strain& BaseElement::strain() const noexcept
{
    return strain_;
}

const Stress& BaseElement::stress() const noexcept
{
    return stress_;
}

bool BaseElement::is_failed() const
{
    return is_failed_;
}


/*****************************************************************************/
/*                              StructElement                                */
/*****************************************************************************/
StructElement::StructElement(int id, Particle& pi, Particle& pj,
    const BaseMaterial& mat, const BaseSection& sect):BaseElement(id)
{
    type_ = "StructElement";
    num_of_particles_ = 2;
    pi_ = &pi;
    pj_ = &pj;
    addParticle(&pi);
    addParticle(&pj);
    fi_ = {0};
    fj_ = {0};
    material_ = &mat;
    section_ = &sect;

    // WARNING: the order can't be changed,
    // or element couldn't be initiated correctly
    calcLength();
    if (length_ <= 0.0)
    {
        throw std::invalid_argument("structural element length must be positive");
    }
    calcMass();
    lt_ = length_;
    la_ = length_;
}

inline void StructElement::calcLength()
{
    const auto& ci = pi_->coordinate();
    const auto& cj = pj_->coordinate();
    length_ = euclideanDistance(ci, cj);
}

inline void StructElement::calcMass()
{
    mass_ = length_ * material_->density() * section_->A();
}

double StructElement::length() const
{
    return length_;
}

double StructElement::calcCriticalTimeStep()
{
    return length_ / material_->calcCe();
}

void StructElement::setParticleForce()
{
    calcElementForce();
    // one particle may be jointed with several element, use addForce()
    pi_->addForce(fi_);
    pj_->addForce(fj_);
}

void StructElement::outputForce(fstream &fout)
{
    fout << id_ << ", ";
    fout << element_force_.fx << ", " << element_force_.sfy << ", ";
    fout << element_force_.sfz << ", " << element_force_.mx << ", ";
    fout << element_force_.my << ", " << element_force_.mz << "\n";
}

void StructElement::outputStress(fstream &fout)
{
    fout << id_ << ", ";
    fout << stress_.sx << ", " << stress_.sy << ", " << stress_.sz << ", ";
    fout << stress_.sxy << ", " << stress_.syz << ", " << stress_.sxz << "\n";
}

void StructElement::outputStrain(fstream &fout)
{
    fout << id_ << ", ";
    fout << strain_.sx << ", " << strain_.sy << ", " << strain_.sz << ", ";
    fout << strain_.sxy << ", " << strain_.syz << ", " << strain_.sxz << "\n";
}

const ElementForce& StructElement::element_force() const noexcept
{
    return element_force_;
}

double StructElement::current_length() const
{
    return lt_;
}


/*****************************************************************************/
/*                               PlaneElement                                */
/*****************************************************************************/
PlaneElement::PlaneElement(int id, vector<Particle*> ps, BaseMaterial* mat,
    BaseSection* sect): BaseElement(id)
{
    type_ = "PlaneElement";
    for (auto iter = ps.begin(); iter != ps.end(); iter++)
    {
        addParticle(*iter);
    }
    material_ = mat;
    section_ = sect;
    // TODO:
}


/*****************************************************************************/
/*                               SolidElement                                */
/*****************************************************************************/
SolidElement::SolidElement(int id, vector<Particle*> ps,
    const BaseMaterial* mat, const BaseSection* sect): BaseElement(id)
{
    type_ = "SolidElement";
    for (auto iter = ps.begin(); iter != ps.end(); iter++)
    {
        addParticle(*iter);
    }
    material_ = mat;
    section_ = sect;
    // TODO:
}


/*****************************************************************************/
/*                                   Link2D                                  */
/*****************************************************************************/
Link2D::Link2D(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
    const BaseSection& sect):StructElement(id, pi, pj, mat, sect)
{
    type_ = "Link2D";
    posi_ = &pi_->current_position();
    posj_ = &pj_->current_position();

    calcOrientVector();
    setParticleDof();
    setParticleMass();
}

void Link2D::setParticleDof()
{
    pi_->activateDof(Dof::Ux);
    pi_->activateDof(Dof::Uy);
    pj_->activateDof(Dof::Ux);
    pj_->activateDof(Dof::Uy);
}

void Link2D::setParticleMass()
{
    StdArray6d mass {mass_/2.0, mass_/2.0, 0, 0, 0, 0};
    pi_->addLumpedMass(mass);
    pj_->addLumpedMass(mass);
}

void Link2D::calcOrientVector()
{
    const auto& ci = pi_->coordinate();
    const auto& cj = pj_->coordinate();
    ex_(0) = (cj[0] - ci[0]) / length_;
    ex_(1) = (cj[1] - ci[1]) / length_;
}

void Link2D::calcElementForce()
{
    calcOrientVector();
    lt_ = euclideanDistance(*posi_, *posj_);
    strain_.sx = (lt_ - la_) / length_;
    stress_.sx = material_->E() * strain_.sx;
    element_force_.fx = stress_.sx * section_->A();

    // fill fi_ and fj_
    fi_.fill(0);
    fj_.fill(0);
    fi_[0] = element_force_.fx * ex_(0);
    fi_[1] = element_force_.fx * ex_(1);
    fj_[0] = -element_force_.fx * ex_(0);
    fj_[1] = -element_force_.fx * ex_(1);
}

void Link2D::outputForce(fstream &fout)
{
    fout << id_ << ", " << element_force_.fx << "\n";
}


/*****************************************************************************/
/*                                   Link2DLD                                */
/*****************************************************************************/
Link2DLD::Link2DLD(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
    const BaseSection& sect): Link2D(id, pi, pj, mat, sect)
{
    type_ = "Link2DLD";
}

void Link2DLD::calcOrientVector()
{
    la_ = lt_;
    lt_ = euclideanDistance(*posi_, *posj_);
    ex_(0) = ((*posj_)[0] - (*posi_)[0]) / lt_;
    ex_(1) = ((*posj_)[1] - (*posi_)[1]) / lt_;
}

void Link2DLD::calcElementForce()
{
    // fill fi_ and fj_
    fi_.fill(0);
    fj_.fill(0);

    // element still working
    if (!is_failed_)
    {
        calcOrientVector();                         // update direction vector
        double delta_strain = (lt_ - la_) / la_;    // increment of strain
        if (material_->isFractured(material_state_.strain + delta_strain))
        {
            // material fractured
            is_failed_ = true;
        }
        else
        {
            // get tagent modulus and update stress and strain status
            material_->update(delta_strain, material_state_);
            strain_.sx = material_state_.strain;
            stress_.sx = material_state_.stress;
            element_force_.fx = material_state_.stress * section_->A();

            // deprecated, same function as previous 3 expressions
            // double delta_stress = eq * delta_strain;
            // double delta_force = delta_stress * section_->A();
            // strain_.sx += delta_strain;
            // stress_.sx += delta_stress;
            // element_force_.fx += delta_force;

            // fill fi_ and fj_
            fi_[0] = element_force_.fx * ex_(0);
            fi_[1] = element_force_.fx * ex_(1);
            fj_[0] = -element_force_.fx * ex_(0);
            fj_[1] = -element_force_.fx * ex_(1);
        }
    }
}


/*****************************************************************************/
/*                                   Link3D                                  */
/*****************************************************************************/
Link3D::Link3D(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
    const BaseSection& sect): StructElement(id, pi, pj, mat, sect)
{
    type_ = "Link3D";
    posi_ = &pi_->current_position();
    posj_ = &pj_->current_position();

    calcOrientVector();
    setParticleDof();
    setParticleMass();
}

void Link3D::setParticleDof()
{
    pi_->activateDof(Dof::Ux);
    pi_->activateDof(Dof::Uy);
    pi_->activateDof(Dof::Uz);
    pj_->activateDof(Dof::Ux);
    pj_->activateDof(Dof::Uy);
    pj_->activateDof(Dof::Uz);
}

void Link3D::setParticleMass()
{
    StdArray6d mass {mass_/2.0, mass_/2.0, mass_/2.0, 0, 0, 0};
    pi_->addLumpedMass(mass);
    pj_->addLumpedMass(mass);
}

void Link3D::calcOrientVector()
{
    const auto& ci = pi_->coordinate();
    const auto& cj = pj_->coordinate();
    ex_(0) = (cj[0] - ci[0]) / length_;
    ex_(1) = (cj[1] - ci[1]) / length_;
    ex_(2) = (cj[2] - ci[2]) / length_;
}

void Link3D::calcElementForce()
{
    calcOrientVector();
    lt_ = euclideanDistance(*posi_, *posj_);
    strain_.sx = (lt_ - la_) / length_;
    stress_.sx = material_->E() * strain_.sx;
    element_force_.fx = stress_.sx * section_->A();

    // fill fi_ and fj_
    fi_.fill(0);
    fj_.fill(0);
    fi_[0] = element_force_.fx * ex_(0);
    fi_[1] = element_force_.fx * ex_(1);
    fi_[2] = element_force_.fx * ex_(2);
    fj_[0] = -element_force_.fx * ex_(0);
    fj_[1] = -element_force_.fx * ex_(1);
    fj_[2] = -element_force_.fx * ex_(2);
}

void Link3D::outputForce(fstream &fout)
{
    fout << id_ << ", " << element_force_.fx << "\n";
}


/*****************************************************************************/
/*                                   Link3DLD                                */
/*****************************************************************************/
Link3DLD::Link3DLD(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
    const BaseSection& sect):Link3D(id, pi, pj, mat, sect)
{
    type_ = "Link3DLD";
}

void Link3DLD::calcOrientVector()
{
    la_ = lt_;
    lt_ = euclideanDistance(*posi_, *posj_);
    ex_(0) = ((*posj_)[0] - (*posi_)[0]) / lt_;
    ex_(1) = ((*posj_)[1] - (*posi_)[1]) / lt_;
    ex_(2) = ((*posj_)[2] - (*posi_)[2]) / lt_;
}

void Link3DLD::calcElementForce()
{
    // fill fi_ and fj_
    fi_.fill(0);
    fj_.fill(0);

    // element still working
    if (!is_failed_)
    {
        calcOrientVector();                         // update direction vector
        double delta_strain = (lt_ - la_) / la_;    // increment of strain
        if (material_->isFractured(material_state_.strain + delta_strain))
        {
            // material fractured
            is_failed_ = true;
        }
        else
        {
            // get tangent modulus and update stress and strain status
            material_->update(delta_strain, material_state_);
            strain_.sx = material_state_.strain;
            stress_.sx = material_state_.stress;
            element_force_.fx = material_state_.stress * section_->A();

            // fill fi_ and fj_
            fi_[0] = element_force_.fx * ex_(0);
            fi_[1] = element_force_.fx * ex_(1);
            fi_[2] = element_force_.fx * ex_(2);
            fj_[0] = -element_force_.fx * ex_(0);
            fj_[1] = -element_force_.fx * ex_(1);
            fj_[2] = -element_force_.fx * ex_(2);
        }
    }
}


/*****************************************************************************/
/*                                   Beam2D                                  */
/*****************************************************************************/
Beam2D::Beam2D(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
    const BaseSection& sect): StructElement(id, pi, pj, mat, sect)
{
    type_ = "Beam2D";
    curr_posi_ = &pi_->current_position();
    curr_posj_ = &pj_->current_position();
    prev_posi_ = &pi_->previous_position();
    prev_posj_ = &pj_->previous_position();
    last_posi_ = *curr_posi_;
    last_posj_ = *curr_posj_;
    elem_rotation_angle_ = 0.0;
    emfi_ = {0, 0, 0, 0, 0, 0};
    emfj_ = {0, 0, 0, 0, 0, 0};

    initiateStiffnessMatrix();
    calcOrientVector();
    setParticleDof();
    setParticleMass();
}

inline void Beam2D::initiateStiffnessMatrix()
{
    // deformation stiffness for [elongation, phi_i, phi_j]; the element
    // stays linear elastic in the co-rotational frame, geometric
    // nonlinearity enters through the frame update in calcElementForce
    const double ea_l = material_->E() * section_->A() / length_;
    const double ei = material_->E() * section_->Izz();
    stiffness_ << ea_l, 0.0, 0.0,
                  0.0, 4.0 * ei / length_, 2.0 * ei / length_,
                  0.0, 2.0 * ei / length_, 4.0 * ei / length_;
    force_.setZero();
}

void Beam2D::setParticleDof()
{
    pi_->activateDof(Dof::Ux);
    pi_->activateDof(Dof::Uy);
    pi_->activateDof(Dof::Rotz);
    pj_->activateDof(Dof::Ux);
    pj_->activateDof(Dof::Uy);
    pj_->activateDof(Dof::Rotz);
}

void Beam2D::setParticleMass()
{
    // half of the bar translates with each particle; the rotary inertia of
    // a half bar of length L/2 about its end is m*L^2/24 (slender rod)
    const double mrot = mass_ * length_ * length_ / 24.0;
    StdArray6d m {mass_ / 2.0, mass_ / 2.0, 0, 0, 0, mrot};
    pi_->addLumpedMass(m);
    pj_->addLumpedMass(m);
}

void Beam2D::calcOrientVector()
{
    const double dx = (*curr_posj_)[0] - (*curr_posi_)[0];
    const double dy = (*curr_posj_)[1] - (*curr_posi_)[1];
    const double l = sqrt(dx * dx + dy * dy);
    exa_ << dx / l, dy / l, 0.0;
    eza_ << 0.0, 0.0, 1.0;
    eya_ = eza_.cross(exa_);
    // expose the current triad through the BaseElement accessors
    ex_ = exa_;
    ey_ = eya_;
    ez_ = eza_;
}

void Beam2D::calcElementForce()
{
    // chord vectors of the last recorded and the current configuration
    const double pdx = last_posj_[0] - last_posi_[0];
    const double pdy = last_posj_[1] - last_posi_[1];
    const double cdx = (*curr_posj_)[0] - (*curr_posi_)[0];
    const double cdy = (*curr_posj_)[1] - (*curr_posi_)[1];

    la_ = sqrt(pdx * pdx + pdy * pdy);
    lt_ = sqrt(cdx * cdx + cdy * cdy);

    // signed rigid-body rotation of the chord over the increment
    const double dtheta_r = atan2(pdx * cdy - pdy * cdx,
                                  pdx * cdx + pdy * cdy);
    elem_rotation_angle_ += dtheta_r;

    // deformation increments: chord elongation and nodal rotations
    // relative to the rigidly rotated chord
    const double du = lt_ - la_;
    const double dphi_i = ((*curr_posi_)[5] - last_posi_[5]) - dtheta_r;
    const double dphi_j = ((*curr_posj_)[5] - last_posj_[5]) - dtheta_r;

    force_ += stiffness_ * Eigen::Vector3d(du, dphi_i, dphi_j);

    const double n = force_(0);        // axial force, tension positive
    const double mi = force_(1);       // end moment at i, ccw positive
    const double mj = force_(2);       // end moment at j, ccw positive
    const double v = (mi + mj) / lt_;  // shear from moment equilibrium

    calcOrientVector();

    // record element forces for outputs
    strain_.sx = (lt_ - length_) / length_;
    stress_.sx = n / section_->A();
    element_force_.fx = n;
    element_force_.sfy = v;
    element_force_.mz = mi;
    emfi_.fx = n;
    emfi_.sfy = v;
    emfi_.mz = mi;
    emfj_.fx = n;
    emfj_.sfy = v;
    emfj_.mz = mj;

    // restoring forces on the particles: minus the nodal forces acting on
    // the element ends [end i: (-N, +V, Mi); end j: (+N, -V, Mj)] rotated
    // into the global frame
    fi_.fill(0);
    fj_.fill(0);
    fi_[0] = n * exa_(0) - v * eya_(0);
    fi_[1] = n * exa_(1) - v * eya_(1);
    fi_[5] = -mi;
    fj_[0] = -n * exa_(0) + v * eya_(0);
    fj_[1] = -n * exa_(1) + v * eya_(1);
    fj_[5] = -mj;

    // remember this configuration for the next increment
    last_posi_ = *curr_posi_;
    last_posj_ = *curr_posj_;
}

double Beam2D::calcCriticalTimeStep()
{
    // Gershgorin bound of M^-1 K over the local beam stiffness gives a
    // safe overestimate of the highest element frequency; dt_cr = 2/w_max
    const double e = material_->E();
    const double a = section_->A();
    const double izz = section_->Izz();
    const double l = length_;
    const double mt = mass_ / 2.0;           // nodal translational mass
    const double mr = mass_ * l * l / 24.0;  // nodal rotary inertia

    const double row_axial = 2.0 * e * a / l;
    const double row_trans = 24.0 * e * izz / (l * l * l)
                           + 12.0 * e * izz / (l * l);
    const double row_rot = 12.0 * e * izz / (l * l) + 6.0 * e * izz / l;

    const double w2 = max(max(row_axial / mt, row_trans / mt), row_rot / mr);
    return 2.0 / sqrt(w2);
}

const ElementForce& Beam2D::end_force_i() const noexcept
{
    return emfi_;
}

const ElementForce& Beam2D::end_force_j() const noexcept
{
    return emfj_;
}

double Beam2D::rotation_angle() const
{
    return elem_rotation_angle_;
}

void Beam2D::outputForce(fstream &fout)
{
    // EID, Ni, Vi, Mi, Nj, Vj, Mj
    fout << id_ << ", ";
    fout << emfi_.fx << ", " << emfi_.sfy << ", " << emfi_.mz << ", ";
    fout << emfj_.fx << ", " << emfj_.sfy << ", " << emfj_.mz << "\n";
}

void Beam2D::outputStress(fstream &fout)
{
    StructElement::outputStress(fout);
}

void Beam2D::outputStrain(fstream &fout)
{
    StructElement::outputStrain(fout);
}


/*****************************************************************************/
/*                                   Beam3D                                  */
/*****************************************************************************/
Beam3D::Beam3D(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
    const BaseSection& sect)
    : Beam3D(id, pi, pj, mat, sect, StdArray3d{0.0, 0.0, 1.0})
{

}

Beam3D::Beam3D(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
    const BaseSection& sect, const StdArray3d ref_node)
    : StructElement(id, pi, pj, mat, sect)
{
    type_ = "Beam3D";
    reference_node_ = ref_node;
    curr_posi_ = &pi_->current_position();
    curr_posj_ = &pj_->current_position();
    prev_posi_ = &pi_->previous_position();
    prev_posj_ = &pj_->previous_position();
    emfi_ = {0, 0, 0, 0, 0, 0};
    emfj_ = {0, 0, 0, 0, 0, 0};

    initiateBeam3D();
    calcOrientVector();
    setParticleDof();
    setParticleMass();
}

inline void Beam3D::initiateBeam3D()
{
    // initial chord direction, read from the same source as
    // Beam2D::calcOrientVector (current_position, which equals coordinate
    // at construction time)
    Eigen::Vector3d Xi((*curr_posi_)[0], (*curr_posi_)[1], (*curr_posi_)[2]);
    Eigen::Vector3d Xj((*curr_posj_)[0], (*curr_posj_)[1], (*curr_posj_)[2]);
    exa_ = (Xj - Xi) / length_;

    // Legacy bepcrd treats the third/reference point as defining the
    // local x-y plane: ez = normalize(ex x r) and ey = ez x ex.  The C++
    // API stores that point/vector directly, so project out any axial
    // component before forming the same right-handed triad.
    Eigen::Vector3d r(reference_node_[0], reference_node_[1],
                       reference_node_[2]);
    Eigen::Vector3d rp = r - r.dot(exa_) * exa_;
    if (rp.norm() < 1e-8)
    {
        r << 0.0, 1.0, 0.0;
        rp = r - r.dot(exa_) * exa_;
    }
    if (rp.norm() < 1e-8)
    {
        r << 1.0, 0.0, 0.0;
        rp = r - r.dot(exa_) * exa_;
    }
    eza_ = exa_.cross(rp).normalized();
    eya_ = eza_.cross(exa_);

    last_posi_ = *curr_posi_;
    last_posj_ = *curr_posj_;

    n_force_ = 0.0;
    t_moment_ = 0.0;
    mz_i_ = 0.0;
    mz_j_ = 0.0;
    mb_i_ = 0.0;
    mb_j_ = 0.0;

    // Seed the stiffness terms from the initial length. calcElementForce()
    // refreshes them from the current/deformed length, matching frm3D.
    const double e = material_->E();
    const double g = material_->G();
    ka_ = e * section_->A() / length_;
    kt_ = g * section_->J() / length_;
    const double eizz_l = e * section_->Izz() / length_;
    const double eiyy_l = e * section_->Iyy() / length_;
    kb1_ << 4.0 * eizz_l, 2.0 * eizz_l,
            2.0 * eizz_l, 4.0 * eizz_l;
    kb2_ << 4.0 * eiyy_l, 2.0 * eiyy_l,
            2.0 * eiyy_l, 4.0 * eiyy_l;
}

void Beam3D::setParticleDof()
{
    pi_->activateDof(Dof::Ux);
    pi_->activateDof(Dof::Uy);
    pi_->activateDof(Dof::Uz);
    pi_->activateDof(Dof::Rotx);
    pi_->activateDof(Dof::Roty);
    pi_->activateDof(Dof::Rotz);
    pj_->activateDof(Dof::Ux);
    pj_->activateDof(Dof::Uy);
    pj_->activateDof(Dof::Uz);
    pj_->activateDof(Dof::Rotx);
    pj_->activateDof(Dof::Roty);
    pj_->activateDof(Dof::Rotz);
}

void Beam3D::setParticleMass()
{
    // xmass.f90 assigns half the element mass to every translational DOF.
    const double half_mass = mass_ / 2.0;
    StdArray6d m {half_mass, half_mass, half_mass, 0.0, 0.0, 0.0};
    pi_->addLumpedMass(m);
    pj_->addLumpedMass(m);

    // The legacy rotational mass is not an isotropic mL^2/24 scalar.  Its
    // principal values are m_half * [J, Iyy, Izz] / A, then projected into
    // global coordinates by the current principal basis.  Store this
    // element contribution so later basis rotations can update each shared
    // particle by matrix delta rather than erasing other elements' inertia.
    Eigen::Matrix3d principal = Eigen::Matrix3d::Zero();
    principal(0, 0) = half_mass * section_->J() / section_->A();
    principal(1, 1) = half_mass * section_->Iyy() / section_->A();
    principal(2, 2) = half_mass * section_->Izz() / section_->A();

    Eigen::Matrix3d basis;
    basis.col(0) = exa_;
    basis.col(1) = eya_;
    basis.col(2) = eza_;
    rotational_mass_contribution_ = basis * principal * basis.transpose();
    pi_->addInertiaMass(rotational_mass_contribution_);
    pj_->addInertiaMass(rotational_mass_contribution_);
}

void Beam3D::calcOrientVector()
{
    // unlike Beam2D this cannot be a memoryless recompute from the current
    // chord alone: the transverse axes carry accumulated twist history, so
    // the actual triad update lives inline in calcElementForce(); this
    // just exposes the current triad through the BaseElement accessors,
    // called once from the constructor and once per step after the triad
    // update (same call sites Beam2D uses its own calcOrientVector from)
    ex_ = exa_;
    ey_ = eya_;
    ez_ = eza_;
}

void Beam3D::calcElementForce()
{
    Eigen::Vector3d Xi((*curr_posi_)[0], (*curr_posi_)[1], (*curr_posi_)[2]);
    Eigen::Vector3d Xj((*curr_posj_)[0], (*curr_posj_)[1], (*curr_posj_)[2]);
    Eigen::Vector3d lastXi(last_posi_[0], last_posi_[1], last_posi_[2]);
    Eigen::Vector3d lastXj(last_posj_[0], last_posj_[1], last_posj_[2]);

    // Legacy frm3D evaluates one displacement increment in the principal
    // triad carried from the previous configuration.  Force evaluation uses
    // that old triad; pcrd is updated only after the new element forces have
    // been assembled.
    la_ = (lastXj - lastXi).norm();
    lt_ = (Xj - Xi).norm();
    const Eigen::Vector3d ex0 = (lastXj - lastXi) / la_;
    const Eigen::Vector3d exd = (Xj - Xi) / lt_;

    Eigen::Vector3d dri((*curr_posi_)[3] - last_posi_[3],
                         (*curr_posi_)[4] - last_posi_[4],
                         (*curr_posi_)[5] - last_posi_[5]);
    Eigen::Vector3d drj((*curr_posj_)[3] - last_posj_[3],
                         (*curr_posj_)[4] - last_posj_[4],
                         (*curr_posj_)[5] - last_posj_[5]);

    // frm3D b3: project nodal rotations into the *old* principal triad.
    const double prx1 = exa_.dot(dri);
    const double pry1 = eya_.dot(dri);
    const double prz1 = eza_.dot(dri);
    const double prx2 = exa_.dot(drj);
    const double pry2 = eya_.dot(drj);
    const double prz2 = eza_.dot(drj);

    // frm3D c: rigid-body bending rotation is the minimal chord rotation.
    // The legacy code uses asin(|ex0 x exd|) and the normalized cross
    // product, so preserve that convention rather than atan2.
    Eigen::Matrix3d Rc = Eigen::Matrix3d::Identity();
    Eigen::Vector3d rb = Eigen::Vector3d::Zero();
    Eigen::Vector3d cross_e = ex0.cross(exd);
    const double cross_norm = cross_e.norm();
    double theta = 0.0;
    // Legacy crossn() resolves chord rotation down to 1e-12 for the
    // deformation split, while bm3D only rotates the carried pcrd basis once
    // |thetan| reaches 1e-10.  Keep those two thresholds distinct: otherwise
    // very small transient bending is incorrectly treated as pure nodal
    // rotation and the first rotational response is lost.
    if (cross_norm > 1e-12)
    {
        const Eigen::Vector3d axis = cross_e / cross_norm;
        theta = asin(std::min(1.0, cross_norm));
        rb = theta * axis;
        if (fabs(theta) >= 1e-10)
            Rc = Eigen::AngleAxisd(theta, axis).toRotationMatrix();
    }

    const double prby = eya_.dot(rb);
    const double prbz = eza_.dot(rb);

    // frm3D d: pure deformations in the old principal directions.
    const double delta = lt_ - la_;
    const double rxd = prx2 - prx1;
    const double ry1d = pry1 - prby;
    const double ry2d = pry2 - prby;
    const double rz1d = prz1 - prbz;
    const double rz2d = prz2 - prbz;

    // frm3D e: all incremental stiffnesses use the current/deformed length.
    const double e = material_->E();
    const double gmat = material_->G();
    ka_ = e * section_->A() / lt_;
    kt_ = gmat * section_->J() / lt_;
    const double eizz_l = e * section_->Izz() / lt_;
    const double eiyy_l = e * section_->Iyy() / lt_;

    n_force_ += ka_ * delta;
    t_moment_ += kt_ * rxd;

    // Legacy bz is bending about local z.  C++ mz has the same sign.
    mz_i_ += (4.0 * rz1d + 2.0 * rz2d) * eizz_l;
    mz_j_ += (2.0 * rz1d + 4.0 * rz2d) * eizz_l;

    // Legacy by is bending about local y.  C++ mb uses the opposite sign so
    // that the restoring-force convention remains Fext + Frestore = 0.
    mb_i_ -= (4.0 * ry1d + 2.0 * ry2d) * eiyy_l;
    mb_j_ -= (2.0 * ry1d + 4.0 * ry2d) * eiyy_l;

    // Assemble restoring nodal loads in the old principal triad.  Negating
    // these gives legacy elf/Fint, which uses Fext - Fint in VFIFE.
    const double v1 = (mz_i_ + mz_j_) / lt_;
    const double v2 = (mb_i_ + mb_j_) / lt_;
    Eigen::Vector3d fi_vec = n_force_ * exa_ - v1 * eya_ - v2 * eza_;
    Eigen::Vector3d fj_vec = -fi_vec;
    Eigen::Vector3d mi_vec =
        -mz_i_ * eza_ + mb_i_ * eya_ + t_moment_ * exa_;
    Eigen::Vector3d mj_vec =
        -mz_j_ * eza_ + mb_j_ * eya_ - t_moment_ * exa_;

    for (int k = 0; k < 3; ++k)
    {
        fi_[k] = fi_vec(k);
        fj_[k] = fj_vec(k);
        fi_[k + 3] = mi_vec(k);
        fj_[k + 3] = mj_vec(k);
    }

    strain_.sx = (lt_ - length_) / length_;
    stress_.sx = n_force_ / section_->A();
    element_force_.fx = n_force_;
    element_force_.sfy = v1;
    element_force_.sfz = v2;
    element_force_.mx = t_moment_;
    element_force_.my = mb_i_;
    element_force_.mz = mz_i_;

    emfi_.fx = n_force_;
    emfi_.sfy = v1;
    emfi_.sfz = v2;
    emfi_.mx = t_moment_;
    emfi_.my = mb_i_;
    emfi_.mz = mz_i_;
    emfj_.fx = n_force_;
    emfj_.sfy = v1;
    emfj_.sfz = v2;
    emfj_.mx = t_moment_;
    emfj_.my = mb_j_;
    emfj_.mz = mz_j_;

    // frm3D updates pcrd only after the force calculation.  First apply the
    // chord/bending rigid-body rotation, then rotate by the node-i principal
    // torsion prx1 about the old chord.  In particular, a pure node-j twist
    // changes torsional deformation but does not rotate the carried triad.
    Eigen::Vector3d new_ex = Rc * exa_;
    Eigen::Vector3d new_ey = Rc * eya_;
    Eigen::Vector3d new_ez = Rc * eza_;
    if (fabs(prx1) >= 1e-10)
    {
        const Eigen::Matrix3d Rt =
            Eigen::AngleAxisd(prx1, ex0).toRotationMatrix();
        new_ex = Rt * new_ex;
        new_ey = Rt * new_ey;
        new_ez = Rt * new_ez;
    }
    exa_ = new_ex;
    eya_ = new_ey;
    eza_ = new_ez;
    calcOrientVector();

    // asmb.f90 reconstructs rotational inertia from the updated pcrd.
    const double half_mass = mass_ / 2.0;
    Eigen::Matrix3d principal = Eigen::Matrix3d::Zero();
    principal(0, 0) = half_mass * section_->J() / section_->A();
    principal(1, 1) = half_mass * section_->Iyy() / section_->A();
    principal(2, 2) = half_mass * section_->Izz() / section_->A();
    Eigen::Matrix3d basis;
    basis.col(0) = exa_;
    basis.col(1) = eya_;
    basis.col(2) = eza_;
    const Eigen::Matrix3d new_rotational_mass =
        basis * principal * basis.transpose();
    const Eigen::Matrix3d rotational_mass_delta =
        new_rotational_mass - rotational_mass_contribution_;
    pi_->addInertiaMass(rotational_mass_delta);
    pj_->addInertiaMass(rotational_mass_delta);
    rotational_mass_contribution_ = new_rotational_mass;

    last_posi_ = *curr_posi_;
    last_posj_ = *curr_posj_;
}

double Beam3D::calcCriticalTimeStep()
{
    // Gershgorin bound of M^-1 K, generalising Beam2D's row structure to
    // axial + torsion + both bending planes (each divided by the matching
    // lumped mass component), dt_cr = 2 / w_max
    const double e = material_->E();
    const double g = material_->G();
    const double a = section_->A();
    const double izz = section_->Izz();
    const double iyy = section_->Iyy();
    const double j = section_->J();
    const double l = length_;
    const double mt = mass_ / 2.0;             // nodal translational mass
    const double mrx = mt * j / a;              // principal rotational mass x
    const double mry = mt * iyy / a;            // principal rotational mass y
    const double mrz = mt * izz / a;            // principal rotational mass z

    const double row_axial = 2.0 * e * a / l;
    const double row_tors = 2.0 * g * j / l;

    // plane 1 (Izz): translation along ey, rotation about ez
    const double row_trans1 = 24.0 * e * izz / (l * l * l)
                             + 12.0 * e * izz / (l * l);
    const double row_rot1 = 12.0 * e * izz / (l * l) + 6.0 * e * izz / l;

    // plane 2 (Iyy): translation along ez, rotation about ey
    const double row_trans2 = 24.0 * e * iyy / (l * l * l)
                             + 12.0 * e * iyy / (l * l);
    const double row_rot2 = 12.0 * e * iyy / (l * l) + 6.0 * e * iyy / l;

    double w2 = row_axial / mt;
    w2 = max(w2, row_tors / mrx);
    w2 = max(w2, row_trans1 / mt);
    w2 = max(w2, row_rot1 / mrz);
    w2 = max(w2, row_trans2 / mt);
    w2 = max(w2, row_rot2 / mry);

    return 2.0 / sqrt(w2);
}

const ElementForce& Beam3D::end_force_i() const noexcept
{
    return emfi_;
}

const ElementForce& Beam3D::end_force_j() const noexcept
{
    return emfj_;
}

void Beam3D::outputForce(fstream &fout)
{
    // EID, Ni, Vyi, Vzi, Ti, Myi, Mzi, Nj, Vyj, Vzj, Tj, Myj, Mzj
    fout << id_ << ", ";
    fout << emfi_.fx << ", " << emfi_.sfy << ", " << emfi_.sfz << ", ";
    fout << emfi_.mx << ", " << emfi_.my << ", " << emfi_.mz << ", ";
    fout << emfj_.fx << ", " << emfj_.sfy << ", " << emfj_.sfz << ", ";
    fout << emfj_.mx << ", " << emfj_.my << ", " << emfj_.mz << "\n";
}

void Beam3D::outputStress(fstream &fout)
{
    StructElement::outputStress(fout);
}

void Beam3D::outputStrain(fstream &fout)
{
    StructElement::outputStrain(fout);
}

} // namespace openvfife
