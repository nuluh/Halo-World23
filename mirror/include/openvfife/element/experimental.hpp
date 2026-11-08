// openVFIFE public API: experimental and unimplemented element declarations
#pragma once

#include <vector>

#include <openvfife/element/frame.hpp>
#include <openvfife/material/uniaxial.hpp>

namespace openvfife {
class PlaneElement: public BaseElement
{
    // TODO:
    private:

    protected:
        void setParticleDof() override = 0;
        void setParticleMass() override = 0;
        void calcOrientVector() override = 0;
        void calcElementForce() override = 0;

    public:
        PlaneElement(int id, std::vector<Particle*> ps, BaseMaterial* mat,
            BaseSection* sect);
        ~PlaneElement() override = default;
        double calcCriticalTimeStep() override = 0;
        void setParticleForce() override = 0;
        void outputForce(std::fstream &fout) override = 0;
        void outputStress(std::fstream &fout) override = 0;
        void outputStrain(std::fstream &fout) override = 0;
};

class SolidElement: public BaseElement
{
    // TODO:
    private:

    protected:
        void setParticleDof() override = 0;
        void setParticleMass() override = 0;
        void calcOrientVector() override = 0;
        void calcElementForce() override = 0;

    public:
        SolidElement(int id, std::vector<Particle*> ps, const BaseMaterial* mat,
            const BaseSection* sect);
        ~SolidElement() override = default;
        double calcCriticalTimeStep() override = 0;
        void setParticleForce() override = 0;
        void outputForce(std::fstream &fout) override = 0;
        void outputStress(std::fstream &fout) override = 0;
        void outputStrain(std::fstream &fout) override = 0;
};

class Link2DLD:public Link2D
{
    private:
        /* parameters for plasitc material
        * ds - delta strain;
        * stress - current stress,
        * strain - current strain,
        * ps - current plastic strain
        * alpha - back stress,
        * q - hardening funciton,
        */
        UniaxialMaterialState material_state_;

    protected:
        void calcOrientVector() override;
        void calcElementForce() override;

    public:
        Link2DLD(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
                const BaseSection& sect);
        ~Link2DLD() override = default;
};

class Link3D: public StructElement
{
    protected:
        const StdArray6d* posi_;
        const StdArray6d* posj_;

        void setParticleDof() override;
        void setParticleMass() override;
        void calcOrientVector() override;
        void calcElementForce() override;

    public:
        Link3D(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
            const BaseSection& sect);
        ~Link3D() override = default;
        void outputForce(std::fstream &fout) override;
};

class Link3DLD:public Link3D
{
    private:
        /* parameters for plasitc material
        * ds - delta strain;
        * stress - current stress,
        * strain - current strain,
        * ps - current plastic strain
        * alpha - back stress,
        * q - hardening funciton,
        */
        UniaxialMaterialState material_state_;

    protected:
        void calcOrientVector() override;
        void calcElementForce() override;

    public:
        Link3DLD(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
               const BaseSection& sect);
        ~Link3DLD() override = default;
};

} // namespace openvfife
