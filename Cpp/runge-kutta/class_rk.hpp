#ifndef CLASS_RK_HPP
#define CLASS_RK_HPP

#include <string>
#include <boost/phoenix/core.hpp>
#include <cmath>
#include <vector>
#include <fstream>
#include <utility>
#include <iostream>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/algorithm/string.hpp>

typedef std::vector<double> state_type;
typedef boost::numeric::ublas::vector<double> vector_type;

// Starting and ending values of lambda
const double x_0 = 1.0;
const double x_c = 0.7;

// Integration step
const double step = -0.0001;

// Astronomical unit
const double au = 150e9;

// Year in seconds
const double year_in_sec = 3600 * 24 * 365;


double sqr(const double &x);

class Bubble {
    private:
        const double gamma, k_rho, n_int, eta;
        const int i, j, k;

    public:
        explicit Bubble(double gamma, double k_rho, double n_int, int i = 0, int j = 0, int k = 0);
        auto CurveValue(const boost::phoenix::placeholders::arg2_type lambda);
        double LambdaApprox();
        std::string Filename();
        /*
        double GradVelShock(double gamma_sa, double k_rho);
        double GradVelCon(double gamma_sa, double k_rho);
        */
        double LambdaShockWind(double lambda_c, double numb_of_years = 1e4, double vel_in = 75000, double r_s = 0);
};

#endif
