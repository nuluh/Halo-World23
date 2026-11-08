// openVFIFE public API: experimental section definitions
#pragma once

#include <openvfife/section/section.hpp>

namespace openvfife {
class Ctube: public BaseSection
{
    /* circular tube section
    *  ri - internal diameter
    *  ro - outer diameter
    *
    *           ^ y
    *           |
    *        ...|...
    *       . . | . .
    *      . .  |  . .
    * ----.-.---|---.-.---------> z
    *      . .  |  . .
    *       . ..|.. .
    *        ...|...
    *           |
    */
    private:
        double ri_;
        double ro_;

    protected:
        void calcArea() override;
        void calcSy() override;
        void calcSz() override;
        void calcCGy() override;
        void calcCGz() override;
        void calcSHy() override;
        void calcSHz() override;
        void calcIyy() override;
        void calcIzz() override;
        void calcIyz() override;

    public:
        Ctube(int id, double ri, double ro);
        ~Ctube() override = default;
        double ri() const;
        double ro() const;
        void info() const override;
};

class AngleSteel: public BaseSection
{
    /* Section Parameters
    * Sy, Sz - calculated in xOy coordinate system
    * CGy, CGZ - presented in xOy coordiante system
    * Iyy, Izz - calculated in xcOyc coordinate system (principle central axis)
    *
    *     y ^ t1
    *       |____
    *       |   |
    *       |   |
    *       |   |    ^ yc
    *  w1   |   |    |
    *       |A1 |    |
    *       |   |    |
    *       |   | ---|-----> zc
    *       |   |    |shape center
    *       |   |______________
    *       |   |     A2      | t2
    *   ____|___|_____________|________> z
    *       |
    *       |        w2
    *
    */
    private:
        double w1_;
        double w2_;
        double t1_;
        double t2_;

    protected:
        void calcArea() override;
        void calcSy() override;
        void calcSz() override;
        void calcCGy() override;
        void calcCGz() override;
        void calcSHy() override;
        void calcSHz() override;
        void calcIyy() override;
        void calcIzz() override;
        void calcIyz() override;

    public:
        AngleSteel(int id, double w1, double w2, double t1, double t2);
        ~AngleSteel() override = default;
        double w1() const;
        double w2() const;
        double t1() const;
        double t2() const;
        void info() const override;
};

} // namespace openvfife
