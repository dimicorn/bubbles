#include <string>
#include <boost/phoenix/core.hpp>
#ifndef CLASS_R_HPP
#define CLASS_R_HPP

class Bubble {
    private:
        double gamma__, k_rho__, n_int__;
        double eta__ = (2 + n_int__) / (5 - k_rho__);
        int i_, j_, k_;

    public:
        explicit Bubble(double gamma, double k_rho, double n_int, int i, int j, int k);
        auto CurveValue(const boost::phoenix::placeholders::arg2_type lambda);
        double LambdaApprox();
        std::string Filename();
        /*
        double GradVel1();
        double GradVel2();
        double Delta(double lambda_c, double lambda_c_approx);
        double R_sw(double vel, double numb, double k);
        */
};

#endif