// openVFIFE public API: system ownership and solver orchestration
#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <openvfife/core/particle.hpp>
#include <openvfife/element/base.hpp>
#include <openvfife/material/uniaxial.hpp>
#include <openvfife/section/section.hpp>

namespace openvfife {
class StructSystem
{
    public:
        using ParticleMap = std::map<int, std::unique_ptr<Particle>>;
        using MaterialMap = std::map<int, std::unique_ptr<BaseMaterial>>;
        using SectionMap = std::map<int, std::unique_ptr<BaseSection>>;
        using ElementMap = std::map<int, std::unique_ptr<BaseElement>>;

    private:
        double h_;
        double alpha_; // mass damping coefficient
        double beta_; // stiffness damping coefficient, default to be 0,

    protected:
        int id_;
        double critical_time_step_;
        std::string workdir_;
        std::string jobname_;
        ParticleMap particles_;
        MaterialMap materials_;
        SectionMap sections_;
        ElementMap elements_;
        std::vector <std::string> elem_types_;
        int elem_max_nodes_;
        std::map <int, DOF> constraints_;
        void constraint();
        void clearOwnedObjects() noexcept;

    public:
        explicit StructSystem(int id);
        ~StructSystem();
        StructSystem(const StructSystem&) = delete;
        StructSystem& operator=(const StructSystem&) = delete;
        StructSystem(StructSystem&&) = delete;
        StructSystem& operator=(StructSystem&&) = delete;
        std::string jobname() const;
        std::string workdir() const;
        double critical_time_step() const;
        double timestep() const;

        [[nodiscard]] const ParticleMap& particles() const noexcept;
        [[nodiscard]] const MaterialMap& materials() const noexcept;
        [[nodiscard]] const SectionMap& sections() const noexcept;
        [[nodiscard]] const ElementMap& elements() const noexcept;
        [[nodiscard]] const std::map<int, DOF>& constraints() const noexcept;
        void info() const;

        void setWorkdir(const std::string &dir);
        void setJobname(const std::string &job);
        double autoTimeStep();
        void setTimeStep(double h);
        void setDampCoeff(double alpha, double beta=0.0);

        void addParticle(std::unique_ptr<Particle> p);
        void addParticles(std::vector<std::unique_ptr<Particle>> particles);

        void addElement(std::unique_ptr<BaseElement> e);
        void addElements(std::vector<std::unique_ptr<BaseElement>> elements);

        void addMaterial(std::unique_ptr<BaseMaterial> mat);
        void addMaterials(std::vector<std::unique_ptr<BaseMaterial>> materials);

        void addSection(std::unique_ptr<BaseSection> sect);
        void addSections(std::vector<std::unique_ptr<BaseSection>> sections);

        void addConstraint(int pid, const DOF &dof);
        [[nodiscard]] bool deleteConstraint(int pid);
        void changeConstraint(int pid, const DOF &dof);
        void clearConstraint();

        void setParticleDisplace(int pid, const StdArray6d &val);
        void addParticleDisplace(int pid, const StdArray6d &val);
        void setAccelerate(double ax, double ay, double az);
        void setExternalForce(int pid, const StdArray6d &val);
        void addExternalForce(int pid, const StdArray6d &val);
        void setInternalForce();
        void clearParticleForce();

        void solve(double h, double zeta=0.2, bool firstStep=false);
        void solve(bool firstStep=false);
        void saveModel(const std::string &path);
        void saveParticleResult(const std::string &path);
        void saveElementResult(const std::string &path);
        void saveSupportReact(const std::string &path);
        void saveResult(const std::string &path);


    // methods for collapse analysis
    protected:
        ElementMap failure_elems_;
        std::map<int, int> particle_connected_elems_;
        ParticleMap free_particles_;

    public:
        void killElement(const int id);
        void autoKillElement();
        void removeFreeParticle();

};

} // namespace openvfife
