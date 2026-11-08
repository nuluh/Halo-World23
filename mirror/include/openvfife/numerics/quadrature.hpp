// openVFIFE public API: spatial quadrature
#pragma once

#include <vector>


namespace openvfife {
class BaseIntegrator
{
    /*
    Atrributes:
        * nodes_ - number of integral nodes
        * x_ - coordinates of integral nodes
        * w_ - the corresponding weights

    Methods:
        integrate1d(const vector<double> &val, double x1, double x2)
        * single integrals
        * val - function values at integral nodes, [f(x1), f(x2), ... f(xn)],
        *       val.size() == nodes
        * x1, x2 - lower and upper limit of the integral

        integrate2d(const vector<vectro<double>> &val, double x1, double x2,
                    double y1, double y2)
        * double integrals
        * val - fuction values at integral nodes, [[f(x1, y1), f(x2, y1),...],
        *                                          [f(x1, y2), f(x2, y2),...],
        *                                          ...
        *                                          [f(x1, yn), f(x2, yn),...]]

        integrate3d(const vector<vector<vector<double>>> &val, double x1,
                    double x2, double y1, double y2, double z1, double z2)
        * triple integrals

    */
    protected:
        int nodes_;
        std::vector<double> x_;
        std::vector<double> w_;

    public:
        BaseIntegrator(int nodes);
        virtual ~BaseIntegrator() = default;
        double x(int i) const;
        double w(int i) const;
        // single integrals
        double integrate1d(const std::vector<double> &val,
                           double x1, double x2) const;
        // double integrals
        double integrate2d(const std::vector<std::vector<double>> &val,
                           double x1, double x2, double y1, double y2) const;
        // triple integrals
        double integrate3d(const std::vector<std::vector<std::vector<double>>>
                           &val, double x1, double x2, double y1, double y2,
                           double z1, double z2) const;
};


class GaussLegendre:public BaseIntegrator
{
    /* A simple version of Gauss-Legendre integration, which support 1~10
    * integartion points at most.
    */
    public:
        GaussLegendre(int nodes);
        ~GaussLegendre() override = default;
};


class GaussLobatto:public BaseIntegrator
{
    /* A simple version of Gauss-Lobatto integration, which support 2~10
    * integartion points at most.
    */
    public:
        GaussLobatto(int nodes);
        ~GaussLobatto() override = default;
};

} // namespace openvfife
