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
        virtual void setParticleDof()=0;
        virtual void setParticleMass()=0;
        virtual void calcOrientVector()=0;
        virtual void calcElementForce()=0;

    public:
        PlaneElement(int id, std::vector<Particle*> ps, BaseMaterial* mat,
            BaseSection* sect);
        ~PlaneElement() override = default;
        virtual double calcCriticalTimeStep()=0;
        virtual void setParticleForce()=0;
        virtual void outputForce(std::fstream &fout)=0;
        virtual void outputStress(std::fstream &fout)=0;
        virtual void outputStrain(std::fstream &fout)=0;
};

class SolidElement: public BaseElement
{
    // TODO:
    private:

    protected:
        virtual void setParticleDof()=0;
        virtual void setParticleMass()=0;
        virtual void calcOrientVector()=0;
        virtual void calcElementForce()=0;

    public:
        SolidElement(int id, std::vector<Particle*> ps, const BaseMaterial* mat,
            const BaseSection* sect);
        ~SolidElement() override = default;
        virtual double calcCriticalTimeStep()=0;
        virtual void setParticleForce()=0;
        virtual void outputForce(std::fstream &fout)=0;
        virtual void outputStress(std::fstream &fout)=0;
        virtual void outputStrain(std::fstream &fout)=0;
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
        virtual void calcOrientVector();
        virtual void calcElementForce();

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

        virtual void setParticleDof();
        virtual void setParticleMass();
        virtual void calcOrientVector();
        virtual void calcElementForce();

    public:
        Link3D(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
            const BaseSection& sect);
        ~Link3D() override = default;
        virtual void outputForce(std::fstream &fout);
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
        virtual void calcOrientVector();
        virtual void calcElementForce();

    public:
        Link3DLD(int id, Particle& pi, Particle& pj, const BaseMaterial& mat,
               const BaseSection& sect);
        ~Link3DLD() override = default;
};

} // namespace openvfife
