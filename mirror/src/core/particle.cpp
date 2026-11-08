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
// Purpose: This file contains the class definition for Particle

// The interface:
//

#include <openvfife/core/particle.hpp>

#include <stdexcept>
#include <string_view>
#include <iostream>
using namespace std;

namespace openvfife {


/*****************************************************************************/
/*                              Helper Functions                             */
/*****************************************************************************/
inline void replaceStdArray(StdArray6d &arr1, const StdArray6d &arr2)
{
    for (int i = 0; i < 6; i++)
    {
        arr1[i] = arr2[i];
    }
};

inline void addStdArray(StdArray6d &arr1, const StdArray6d &arr2)
{
    for (int i = 0; i < 6; i++)
    {
        arr1[i] += arr2[i];
    }
};


inline void validateLumpedMass(const StdArray6d &mass)
{
    for (const double value : mass)
    {
        if (value < 0.0)
        {
            throw std::invalid_argument("lumped mass must be non-negative");
        }
    }
}

inline void validateTimeStep(double h)
{
    if (!(h > 0.0))
    {
        throw std::invalid_argument("time step must be positive");
    }
}

inline void validateDamping(double zeta)
{
    if (zeta < 0.0)
    {
        throw std::invalid_argument("damping coefficient must be non-negative");
    }
}

inline void validateMassConfiguration(const Mass &mass,
    bool translate_mass_matrix_on, bool rotation_mass_matrix_on)
{
    if (translate_mass_matrix_on)
    {
        throw std::logic_error("translation mass-matrix integration is not implemented");
    }
    if (rotation_mass_matrix_on &&
        !mass.rotation_mass_matrix.fullPivLu().isInvertible())
    {
        throw std::domain_error("rotation mass matrix must be invertible");
    }
}

constexpr std::array<std::string_view, 6> dof_names{
    "Ux", "Uy", "Uz", "Rotx", "Roty", "Rotz"
};

inline std::size_t dofIndex(Dof dof)
{
    const auto index = static_cast<std::size_t>(dof);
    if (index >= dof_names.size())
    {
        throw std::invalid_argument("unknown DOF");
    }
    return index;
}

/*****************************************************************************/
/*                                  Particle                                 */
/*****************************************************************************/
inline void Particle::initiate()
{
    current_position_ = coordinate_;
    previous_position_ = coordinate_;
}

Particle::Particle(int id)
    : id_(id)
{
    initiate();
}

Particle::Particle(int id, double x, double y)
    : Particle(id, x, y, 0.0)
{
}

Particle::Particle(int id, double x, double y, double z)
    : id_(id), coordinate_{x, y, z, 0.0, 0.0, 0.0}
{
    initiate();
}

void Particle::info() const
{
    cout << "Particle ID: " << id_ << endl;
    cout << "coordinate: (";
    for (int i = 0; i < 6; i++)
    {
        cout << coordinate_[i] << ", ";
    }
    cout << ")" << endl;

    cout << "dofs: ";
    for (std::size_t i = 0; i < dof_index_.size(); ++i)
    {
        if (dof_index_[i])
        {
            cout << dof_names[i] << ", ";
        }
    }
    cout << endl;

    cout << "mass:\n" << "lumped_mass: ";
    for (int i = 0; i < 6; i++)
    {
        cout << mass_.lumped_mass[i] << ", ";
    }
    cout << endl;
    cout << "consistent translate mass:" << endl;
    cout << mass_.translate_mass_matrix << endl;

    cout << "consistent rotation mass:" << endl;
    cout << mass_.rotation_mass_matrix << endl;

    cout << "force: ";
    for (int i = 0; i < 6; i++)
    {
        cout << force_[i] << ", ";
    }
    cout << "\n" << endl;
}

int Particle::id() const
{
    return id_;
}

const StdArray6d& Particle::coordinate() const noexcept
{
    return coordinate_;
}
const Mass& Particle::mass() const noexcept
{
    return mass_;
}

const StdArray6d& Particle::force() const noexcept
{
    return force_;
}

const StdArray6d& Particle::displace() const noexcept
{
    return displace_;
}


const StdArray6d& Particle::velocity() const noexcept
{
    return velocity_;
}

const StdArray6d& Particle::accelerate() const noexcept
{
    return accelerate_;
}


const StdArray6d& Particle::previous_displace() const noexcept
{
    return previous_displace_;
}

const StdArray6d& Particle::previous_velocity() const noexcept
{
    return previous_velocity_;
}

const StdArray6d& Particle::previous_accelerate() const noexcept
{
    return previous_accelerate_;
}

const StdArray6d& Particle::current_position() const noexcept
{
    return current_position_;
}

const StdArray6d& Particle::previous_position() const noexcept
{
    return previous_position_;
}

const StdArray6b& Particle::dof_index() const noexcept
{
    return dof_index_;
}

void Particle::setCoordinate(const StdArray6d &val)
{
    replaceStdArray(coordinate_, val);
}

void Particle::activateDof(Dof dof)
{
    dof_index_[dofIndex(dof)] = true;
}

void Particle::deactivateDof(Dof dof)
{
    dof_index_[dofIndex(dof)] = false;
}

void Particle::constraintDof(const DOF &dof)
{
    for (int i = 0; i < 6; i++)
    {
        if ((dof_index_[i]) && dof.key[i])
        {
            displace_[i] = dof.val[i];
            previous_displace_[i] = dof.val[i];
            velocity_[i] = 0;
            previous_velocity_[i] = 0;
            accelerate_[i] = 0;
            previous_accelerate_[i] = 0;
        }
    }
    updatePosition();
}

void Particle::setLumpedMass(const StdArray6d &val)
{
    validateLumpedMass(val);
    replaceStdArray(mass_.lumped_mass, val);
}

void Particle::addLumpedMass(const StdArray6d &val)
{
    StdArray6d candidate = mass_.lumped_mass;
    addStdArray(candidate, val);
    validateLumpedMass(candidate);
    mass_.lumped_mass = candidate;
}

void Particle::clearLumpedMass()
{
    for (int i = 0; i < 6; i++)
    {
        mass_.lumped_mass[i] = 0;
    }
}

void Particle::setTranslateMass(const Eigen::Matrix3d &mat)
{
    translate_mass_matrix_on_ = true;
    mass_.translate_mass_matrix << mat;
}

void Particle::addTranslateMass(const Eigen::Matrix3d &mat)
{
    translate_mass_matrix_on_ = true;
    mass_.translate_mass_matrix += mat;
}

void Particle::clearTranslateMass()
{
    translate_mass_matrix_on_ = false;
    mass_.translate_mass_matrix.setZero();
}

void Particle::setInertiaMass(const Eigen::Matrix3d &mat)
{
    rotation_mass_matrix_on_ = true;
    mass_.rotation_mass_matrix << mat;
}

void Particle::addInertiaMass(const Eigen::Matrix3d &mat)
{
    rotation_mass_matrix_on_ = true;
    mass_.rotation_mass_matrix += mat;
}

void Particle::clearInertiaMass()
{
    rotation_mass_matrix_on_ = false;
    mass_.rotation_mass_matrix.setZero();
}

void Particle::setForce(const StdArray6d &val)
{
    replaceStdArray(force_, val);
}

void Particle::addForce(const StdArray6d &val)
{
    addStdArray(force_, val);
}

void Particle::clearForce()
{
    for (int i = 0; i < 6; i++)
    {
        force_[i] = 0;
    }
}

void Particle::setDisplace(const StdArray6d &val)
{
    replaceStdArray(displace_, val);
}

void Particle::addDisplace(const StdArray6d &val)
{
    addStdArray(displace_, val);
}

void Particle::clearDisplace()
{
    for (int i = 0; i < 6; i++)
    {
        displace_[i] = 0;
    }
}

void Particle::setVelocity(const StdArray6d &val)
{
    replaceStdArray(velocity_, val);
}

void Particle::addVelocity(const StdArray6d &val)
{
    addStdArray(velocity_, val);
}

void Particle::clearVelocity()
{
    for (int i = 0; i < 6; i++)
    {
        velocity_[i] = 0;
    }
}

void Particle::setAccelerate(const StdArray6d &val)
{
    replaceStdArray(accelerate_, val);
}

void Particle::addAccelerate(const StdArray6d &val)
{
    addStdArray(accelerate_, val);
}

void Particle::clearAccelerate()
{
    for (int i = 0; i < 6; i++)
    {
        accelerate_[i] = 0;
    }
}

void Particle::setPreviousDisplace(const StdArray6d &val)
{
    replaceStdArray(previous_displace_, val);
}

void Particle::addPreviousDisplace(const StdArray6d &val)
{
    addStdArray(previous_displace_, val);
}

void Particle::clearPreviousDisplace()
{
    for (int i = 0; i < 6; i++)
    {
        previous_displace_[i] = 0;
    }
}

void Particle::setPreviousVelocity(const StdArray6d &val)
{
    replaceStdArray(previous_velocity_, val);
}

void Particle::addPreviousVelocity(const StdArray6d &val)
{
    addStdArray(previous_velocity_, val);
}

void Particle::clearPreviousVelocity()
{
    for (int i = 0; i < 6; i++)
    {
        if (dof_index_[i])
        {
            previous_velocity_[i] = 0;
        }
    }
}

void Particle::setPreviousAccelerate(const StdArray6d &val)
{
    replaceStdArray(previous_accelerate_, val);
}

void Particle::addPreviousAccelerate(const StdArray6d &val)
{
    addStdArray(previous_accelerate_, val);
}

void Particle::clearPreviousAccelerate()
{
    for (int i = 0; i < 6; i++)
    {
        previous_accelerate_[i] = 0;
    }
}

void Particle::updatePosition()
{
    for (int i = 0; i < 6; i++)
    {
        current_position_[i] = coordinate_[i] + displace_[i];
        previous_position_[i] = coordinate_[i] + previous_displace_[i];
    }
}

void Particle::clearAll()
{
    clearDisplace();
    clearVelocity();
    clearAccelerate();
    clearPreviousDisplace();
    clearPreviousVelocity();
    clearPreviousAccelerate();
    clearForce();
    // reset position
    for (int i = 0; i < 6; i++)
    {
        current_position_[i] = coordinate_[i];
        previous_position_[i] = coordinate_[i];
    }
}

void Particle::solveFirstStep(double h, double zeta)
{
    validateTimeStep(h);
    validateDamping(zeta);
    validateMassConfiguration(mass_, translate_mass_matrix_on_,
                              rotation_mass_matrix_on_);
    if ((!translate_mass_matrix_on_) && (!rotation_mass_matrix_on_))
    {
        for (int i = 0; i < 6; i++)
        {
            if (dof_index_[i])
            {
                accelerate_[i] = force_[i] / mass_.lumped_mass[i];
                previous_displace_[i] = accelerate_[i] * h * h / 2.0 -
                                    h * (1 + 0.5 * zeta * h) * velocity_[i] +
                                    displace_[i];
            }
        }
    }
    else if (rotation_mass_matrix_on_)
    {
        for (int i = 0; i < 3; i++)
        {
            if (dof_index_[i])
            {
                accelerate_[i] = force_[i] / mass_.lumped_mass[i];
                previous_displace_[i] = accelerate_[i] * h * h / 2.0 -
                                    h * (1 + 0.5 * zeta * h) * velocity_[i] +
                                    displace_[i];
            }
        }

        Eigen::Vector3d moment, tmp;
        moment << force_[3], force_[4], force_[5];
        tmp = mass_.rotation_mass_matrix.inverse() * moment;

        for (int i = 3; i < 6; i++)
        {
            accelerate_[i] = tmp(i-3);
            previous_displace_[i] = accelerate_[i] * h * h / 2.0 -
                                h * (1 + 0.5 * zeta * h) * velocity_[i] +
                                displace_[i];

        }
    }

    updatePosition();
}

void Particle::solve(double h, double zeta)
{
    validateTimeStep(h);
    validateDamping(zeta);
    validateMassConfiguration(mass_, translate_mass_matrix_on_,
                              rotation_mass_matrix_on_);
    const double c1 = 1.0 + zeta * h / 2.0;
    const double c2 = 1.0 - zeta * h / 2.0;
    const double h2 = h * h;

    StdArray6d disp = {};
    StdArray6d velo = {};
    StdArray6d acc = {};

    if ((!translate_mass_matrix_on_) && (!rotation_mass_matrix_on_))
    {
        for (int i = 0; i < 6; i++)
        {
            if (dof_index_[i])
            {
                disp[i] = (2.0 * displace_[i] - c2 * previous_displace_[i] +
                          h2 * force_[i] / mass_.lumped_mass[i]) / c1;
                velo[i] = (disp[i] - previous_displace_[i]) / h / 2.0;
                acc[i] = (disp[i] - 2.0 * displace_[i] +
                          previous_displace_[i]) / h2;
            }
        }
    }
    else if (rotation_mass_matrix_on_)
    {
        for (int i = 0; i < 3; i++)
        {
            if (dof_index_[i])
            {
                disp[i] = (2.0 * displace_[i] - c2 * previous_displace_[i] +
                          h2 * force_[i] / mass_.lumped_mass[i]) / c1;
                velo[i] = (disp[i] - previous_displace_[i]) / h / 2.0;
                acc[i] = (disp[i] - 2.0 * displace_[i] +
                          previous_displace_[i]) / h2;
            }
        }

        Eigen::Vector3d moment, tmp;
        moment << force_[3], force_[4], force_[5];
        tmp = mass_.rotation_mass_matrix.inverse() * moment;

        for (int i = 3; i < 6; i++)
        {
            disp[i] = (2.0 * displace_[i] - c2 * previous_displace_[i] +
                       h2 * tmp(i-3)) / c1;
            velo[i] = (disp[i] - previous_displace_[i]) / h / 2.0;

            acc[i] = (disp[i] - 2.0 * displace_[i] +
                      previous_displace_[i]) / h2;
        }
    }

    setPreviousDisplace(displace_);
    setPreviousVelocity(velocity_);
    setPreviousAccelerate(accelerate_);
    setDisplace(disp);
    setVelocity(velo);
    setAccelerate(acc);
    updatePosition();
}

void Particle::outputForce(fstream &fout) const
{
    fout << id_ << ", ";
    for (int i = 0; i < 5; i++)
    {
        fout << force_[i] << ", ";
    }
    fout << force_[5] << endl;
}

void Particle::outputReactionForce(fstream &fout) const
{
    fout << id_ << ", ";
    for (int i = 0; i < 5; i++)
    {
        fout << -force_[i] << ", ";
    }
    fout << -force_[5] << endl;
}

void Particle::outputDisplace(fstream &fout) const
{
    fout << id_ << ", ";
    for (int i = 0; i < 5; i++)
    {
        fout << displace_[i] << ", ";
    }
    fout << displace_[5] << endl;
}

void Particle::outputVelocity(fstream &fout) const
{
    fout << id_ << ", ";
    for (int i = 0; i < 5; i++)
    {
        fout << velocity_[i] << ", ";
    }
    fout << velocity_[5] << endl;
}

void Particle::outputAccelarate(fstream &fout) const
{
    fout << id_ << ", ";
    for (int i = 0; i < 5; i++)
    {
        fout << accelerate_[i] << ", ";
    }
    fout << accelerate_[5] << endl;
}

} // namespace openvfife
