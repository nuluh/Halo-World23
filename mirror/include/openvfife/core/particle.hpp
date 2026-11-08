// openVFIFE public API: particle and DOF state
#pragma once

#include <cstdint>
#include <fstream>

#include "Eigen/Dense"
#include <openvfife/core/types.hpp>

namespace openvfife {
enum class Dof : std::uint8_t
{
    Ux = 0,
    Uy,
    Uz,
    Rotx,
    Roty,
    Rotz,
};

struct DOF
{
    StdArray6b key{};
    StdArray6d val{};
};

struct Mass
{
    StdArray6d lumped_mass{};
    Eigen::Matrix3d translate_mass_matrix{Eigen::Matrix3d::Zero()}; // mass matrix
    Eigen::Matrix3d rotation_mass_matrix{Eigen::Matrix3d::Zero()};  // inertia mass
};

class Particle
{
    /*
    * control equation of particle, Newton's second law
    * $ m \ddot{x} = f $
    * so, particle should store mass, motion and force information.
    * To solve the equation, the second order centeral difference partial
    * method is adopted, which is also implemented in Particle class, see
    * solve() method.
    *
    *Particle has 6 dofs ['Ux', 'Uy', 'Uz', 'Rotx', 'Roty', 'Rotz'];
    *hence mass_, force_ and displace_ etc. have 6 elements;
    *for mass_ [0, 1, 2, 3, 4, 5] --> [Mx, My, Mz, Imx, Imy, Imz];
    *for Im_ --> [[Imxx, Imxy, Imxz],
    *             [Imyx, Imyy, Imyz],
    *             [Imzx, Imzy, Imzz]]
    *for force_ [0, 1, 2, 3, 4, 5] --> [Fx, Fy, Fz, Mx, My, Mz];
    *for dispalce_ etc. [0, 1, 2, 3, 4, 5] --> [Ux, Uy, Uz, Rotx, Roty, Rotz];
    *
    *dofKey_ stores the dofs of one Particle, (Ux, Uy, Uz, Rotx, Roty, Rotz);
    *dofIndex_ stores the status of each dof, for example:
    *    [true, true, false, false, false, false] means Ux and Uy are activated
    *    however the rest of dofs are deactivated;
    */
    private:
        int id_{};
        StdArray6d coordinate_{};
        Mass mass_{};
        StdArray6d force_{};
        StdArray6d displace_{};
        StdArray6d velocity_{};
        StdArray6d accelerate_{};
        StdArray6d previous_displace_{};
        StdArray6d previous_velocity_{};
        StdArray6d previous_accelerate_{};
        StdArray6d current_position_{};
        StdArray6d previous_position_{};
        StdArray6b dof_index_{};
        bool translate_mass_matrix_on_{false};
        bool rotation_mass_matrix_on_{false};

    private:
        inline void initiate();

    public:
        explicit Particle(int id);
        Particle(int id, double x, double y);
        Particle(int id, double x, double y, double z);
        void info() const;

        int id() const;
        [[nodiscard]] const StdArray6d& coordinate() const noexcept;
        [[nodiscard]] const Mass& mass() const noexcept;
        [[nodiscard]] const StdArray6d& force() const noexcept;

        [[nodiscard]] const StdArray6d& displace() const noexcept;
        [[nodiscard]] const StdArray6d& velocity() const noexcept;
        [[nodiscard]] const StdArray6d& accelerate() const noexcept;

        [[nodiscard]] const StdArray6d& previous_displace() const noexcept;
        [[nodiscard]] const StdArray6d& previous_velocity() const noexcept;
        [[nodiscard]] const StdArray6d& previous_accelerate() const noexcept;

        [[nodiscard]] const StdArray6d& current_position() const noexcept;
        [[nodiscard]] const StdArray6d& previous_position() const noexcept;
        [[nodiscard]] const StdArray6b& dof_index() const noexcept;

        void setCoordinate(const StdArray6d &val);
        void activateDof(Dof dof);
        void deactivateDof(Dof dof);
        void constraintDof(const DOF &dof);

        void setLumpedMass(const StdArray6d &val);
        void addLumpedMass(const StdArray6d &val);
        void clearLumpedMass();

        void setTranslateMass(const Eigen::Matrix3d &mat);
        void addTranslateMass(const Eigen::Matrix3d &mat);
        void clearTranslateMass();

        void setInertiaMass(const Eigen::Matrix3d &mat);
        void addInertiaMass(const Eigen::Matrix3d &mat);
        void clearInertiaMass();

        void setForce(const StdArray6d &val);
        void addForce(const StdArray6d &val);
        void clearForce();

        void setDisplace(const StdArray6d &val);
        void addDisplace(const StdArray6d &val);
        void clearDisplace();

        void setVelocity(const StdArray6d &val);
        void addVelocity(const StdArray6d &val);
        void clearVelocity();

        void setAccelerate(const StdArray6d &val);
        void addAccelerate(const StdArray6d &val);
        void clearAccelerate();

        void setPreviousDisplace(const StdArray6d &val);
        void addPreviousDisplace(const StdArray6d &val);
        void clearPreviousDisplace();

        void setPreviousVelocity(const StdArray6d &val);
        void addPreviousVelocity(const StdArray6d &val);
        void clearPreviousVelocity();

        void setPreviousAccelerate(const StdArray6d &val);
        void addPreviousAccelerate(const StdArray6d &val);
        void clearPreviousAccelerate();

        void updatePosition();
        void clearAll();
        void solveFirstStep(double h, double zeta);
        void solve(double h, double zeta);

        void outputForce(std::fstream &fout) const;
        void outputReactionForce(std::fstream &fout) const;
        void outputDisplace(std::fstream &fout) const;
        void outputVelocity(std::fstream &fout) const;
        void outputAccelarate(std::fstream &fout) const;
};

} // namespace openvfife
