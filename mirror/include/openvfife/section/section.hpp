// openVFIFE public API: supported section definitions
#pragma once

#include <string>

namespace openvfife {
class BaseSection
{
    /* Abstract Base class for Section Module; all class inherited from
    * BaseSection should implement a series method to calculation section
    * characteritics, see all the pure virtual functions.
    *
    */
    protected:
        std::string type_;          // section type name
        int id_;                    // section id
        double A_;                  // area
        double Sy_, Sz_;            // moment of area
        double CGy_, CGz_;          // shape center
        double SHy_, SHz_;          // shear center
        double Iyy_, Izz_, Iyz_;    // inertial moments and product
        double J_;                  // torsion constant (St. Venant)

    protected:
        virtual void calcArea() = 0;
        virtual void calcSy() = 0;
        virtual void calcSz() = 0;
        virtual void calcCGy() = 0;
        virtual void calcCGz() = 0;
        virtual void calcSHy() = 0;
        virtual void calcSHz() = 0;
        virtual void calcIyy() = 0;
        virtual void calcIzz() = 0;
        virtual void calcIyz() = 0;
        void initiate();

    public:
        BaseSection(int id);
        virtual ~BaseSection() = default;
        virtual void info() const;
        std::string type() const;
        int id() const;
        double A() const;
        double Sy() const;
        double Sz() const;
        double CGy() const;
        double CGz() const;
        double SHy() const;
        double SHz() const;
        double Iyy() const;
        double Izz() const;
        double Iyz() const;
        double J() const;
};

struct CustomSectionParas
{
    /* struct for store section characteristics, see CustomSection class */
    double A{};
    double Sy{}, Sz{};
    double SHy{}, SHz{};
    double CGy{}, CGz{};
    double Iyy{}, Izz{}, Iyz{};
    double J{};  // torsion constant; if <= 0 the CustomSection constructor
                 // falls back to the polar approximation J = Iyy + Izz
};

class CustomSection: public BaseSection
{
    /* class for user defined section, all the section characteristics should
    * be inputed by users by using CustomSectionParas struct.
    * Parameter lists:
    *   A - area
    *   Sy, Sz - Moment of area of y and z axis, respectively.
    *   SHy, SHz - Shear center
    *   CGy, CGz - shape center
    *   Iyy, Izz - moment of inertia
    *   Iyz - inertial product
    */
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
        CustomSection(int id, const CustomSectionParas &paras);
        ~CustomSection() override = default;
        void info() const override;
};

class Rectangle: public BaseSection
{
    /* Seciton Parameters
    * xOy coordinate system is same as xcOyc;
    *
    *
    *              ^ y
    *        ______|______
    *        |     |     |
    *        |     |     |
    *        |     |     |  height
    *    ----|-----|-----|-------> z
    *        |     |     |
    *        |     |     |
    *        |_____|_____|
    *              |
    *           width
    *
    */
    protected:
        double width_;
        double height_;

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
        Rectangle(int id, double width, double height);
        ~Rectangle() override = default;
        double width() const;
        double height() const;
        void info() const override;
};

} // namespace openvfife
