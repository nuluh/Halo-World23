// openVFIFE public API: experimental uniaxial materials
#pragma once

#include <openvfife/material/uniaxial.hpp>

namespace openvfife {
class UniIdeal: public UniaxialPlastic
{
    protected:
        double Et_;  // Tangent Modulus, Et_ == 0

    public:
        UniIdeal(int id, double E, double yield_stress);
        ~UniIdeal() override = default;
        // virtual void info() const;
        double update(double delta_strain,
                      UniaxialMaterialState &state) const override;
};

class UniPower: public UniaxialPlastic
{
    // TODO: can't work! The Eq is not right(need to record last yield strain)
    /* stress = E * strain (not yield), stress = k * epsilon ** n (yield)*/
    /*
    * Equation: $\sigma = E \epsilon for \sigma <= {\sigma}_0$
    *           $\sigma = k {\epsilon}^n for \sigma > {\sigma}_0$
    *           ${\sigma}_0 = k({\sigma}_{0} / E) ^ n$
    */
    private:
        double k_;
        double n_;

    protected:

    public:
        UniPower(int id, double E, double n, double yield_stress, double m);
        ~UniPower() override = default;
        virtual void info() const;
        double update(double delta_strain,
                      UniaxialMaterialState &state) const override;
};

} // namespace openvfife
