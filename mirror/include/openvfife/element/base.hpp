// openVFIFE public API: common element contracts
#pragma once

#include <fstream>
#include <map>
#include <string>

#include "Eigen/Dense"
#include <openvfife/core/particle.hpp>
#include <openvfife/material/uniaxial.hpp>
#include <openvfife/section/section.hpp>

namespace openvfife {
struct Strain
{
    /* organize strain at one point, for outputing */
    double sx{}, sy{}, sz{};
    double sxy{}, syz{}, sxz{};
};

struct Stress
{
    /* organize stress at one point, for outputing */
    double sx{}, sy{}, sz{};
    double sxy{}, syz{}, sxz{};
};

struct ElementForce
{
    /* organize force of rod element, for outputing */
    double fx{};
    double sfy{}, sfz{};
    double mx{}, my{}, mz{};
};

/* Sorry, I have no time to write a detailed comments for every element class.
* In fact, I am not very satisfied with the current framework. I will redesign
* openVFIFE after graduating. When the time comes, a detailed documentation
* will be released with new version of openVFIFE. Sorry again, there will be
* no comments of element classes.
*
*/

class BaseElement
{
    protected:
        int id_{};
        std::string type_{"BaseElement"};
        int num_of_particles_{1};
        std::map<int, Particle*> particles_;
        const BaseMaterial* material_{nullptr};
        const BaseSection* section_{nullptr};
        double mass_{};
        Strain strain_{};
        Stress stress_{};
        Eigen::Vector3d ex_{Eigen::Vector3d::Zero()};
        Eigen::Vector3d ey_{Eigen::Vector3d::Zero()};
        Eigen::Vector3d ez_{Eigen::Vector3d::Zero()};
        bool is_failed_{false};

    protected:
        void addParticle(Particle* p);
        virtual void setParticleDof()=0;
        virtual void setParticleMass()=0;
        virtual void calcOrientVector()=0;
        virtual void calcElementForce()=0;

    public:
        explicit BaseElement(int id);
        virtual ~BaseElement() = default;

        virtual void info() const;
        int id() const;
        std::string type() const;
        int num_of_particles() const;
        [[nodiscard]] const std::map<int, Particle*>& particles() const noexcept;
        [[nodiscard]] const BaseMaterial& material() const noexcept;
        [[nodiscard]] const BaseSection& section() const noexcept;
        double mass() const;
        [[nodiscard]] const Eigen::Vector3d& ex() const noexcept;
        [[nodiscard]] const Eigen::Vector3d& ey() const noexcept;
        [[nodiscard]] const Eigen::Vector3d& ez() const noexcept;
        [[nodiscard]] const Strain& strain() const noexcept;
        [[nodiscard]] const Stress& stress() const noexcept;
        bool is_failed() const;
        virtual double calcCriticalTimeStep()=0;
        virtual void setParticleForce()=0;
        virtual void outputForce(std::fstream &fout)=0;
        virtual void outputStress(std::fstream &fout)=0;
        virtual void outputStrain(std::fstream &fout)=0;
};

class StructElement: public BaseElement
{
    private:
        inline void calcLength();
        inline void calcMass();

    protected:
        Particle* pi_{nullptr};
        Particle* pj_{nullptr};
        StdArray6d fi_{};
        StdArray6d fj_{};
        double length_{};
        double lt_{}, la_{};
        ElementForce element_force_{};

        void setParticleDof() override = 0;
        void setParticleMass() override = 0;
        void calcOrientVector() override = 0;
        void calcElementForce() override = 0;

    public:
        StructElement(int id, Particle& pi, Particle& pj,
            const BaseMaterial& mat, const BaseSection& sect);
        ~StructElement() override = default;
        double length() const;
        double current_length() const;
        [[nodiscard]] const ElementForce& element_force() const noexcept;
        double calcCriticalTimeStep() override;
        void setParticleForce() override;
        // virtual void kill()=0;
        void outputForce(std::fstream &fout) override;
        void outputStress(std::fstream &fout) override;
        void outputStrain(std::fstream &fout) override;
};

} // namespace openvfife
