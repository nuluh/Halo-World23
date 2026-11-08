// openVFIFE public API: supported uniaxial materials
#pragma once

#include <string>


namespace openvfife {
struct UniaxialMaterialState
{
    double stress{};
    double strain{};
    double plastic_strain{};
    double back_stress{};
    double hardening{};
};

class BaseMaterial
{
    // when E, G and nu are all defined; the G = E / 2 / (1 + nu);
    protected:
        int id_;
        std::string type_;
        double E_;
        double G_;
        double nu_;
        double density_;
        double ce_;
        double damp_ratio_;
        double damp_coef_;
        double mu_;
        double yield_stress_;
        double limited_tensile_strain_;

    public:
        BaseMaterial(int id);
        virtual ~BaseMaterial() = default;

        int id() const;
        std::string type() const;
        double E() const;
        double G() const;
        double nu() const;
        double density() const;
        double ce() const;
        double damp_ratio() const;
        double damp_coef() const;
        double mu() const;
        double yield_stress() const;
        double limited_tensile_strain() const;
        virtual bool isFractured(double strain) const;
        virtual double update(double delta_strain,
                              UniaxialMaterialState &state) const = 0;

        void setE(double val);
        void setG(double val);
        void setNu(double val);
        void setDensity(double val);
        void setDampRatio(double val);
        void setDampCoef(double val);
        void setMu(double val);
        void setLimitedTensileStrain(double tensile);
        virtual void info() const;

        double calcCe() const;
};

class LinearElastic: public BaseMaterial
{
    public:
        LinearElastic(int id);
        ~LinearElastic() override = default;
        // will change stress and strain
        double update(double delta_strain,
                      UniaxialMaterialState &state) const override;
};

class UniaxialPlastic: public BaseMaterial
{
    // pure virtual class
    protected:
        double m_;  // hardening, 0 - kinematic, 1 - isotropic, 0~1 - mixed

    public:
        UniaxialPlastic(int id);
        ~UniaxialPlastic() override = default;
        void setYieldStress(double val);
        virtual void info() const;
};

class UniBilinear: public UniaxialPlastic
{
    protected:
        double Et_;  // Tagent Modulus

    public:
        UniBilinear(int id, double E, double Et, double yield_stress, double m);
        ~UniBilinear() override = default;
        double Et() const;
        virtual void info() const;
        double update(double delta_strain,
                      UniaxialMaterialState &state) const override;
};

} // namespace openvfife
