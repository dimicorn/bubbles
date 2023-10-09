#pragma once

#include <string>
#include <boost/phoenix/core.hpp>
#include <cmath>
#include <vector>
#include <fstream>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <mpi.h>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/algorithm/string.hpp>

typedef std::vector<double> state_type;
typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

const std::string unknown_arg = "unknown argument: ";
const std::string too_many_args = "too many arguments";

// Starting and ending values of lambda
const double x_0 = 1.0;
const double x_c = 0.7;

// Integration step
const double step = -0.0001;

// Astronomical unit
const double au = 150e9;

// Year in seconds
const double year_in_sec = 3600 * 24 * 365;

// Precision
const double eps = 1.0e-8;

double sqr(const double &x);

enum class Methods {
	ROSENBROCK,
	RUNGE_KUTTA
};

class Bubble {
    private:
        const double gamma, k_rho, n_int, eta;
        const int i, j, k;
		const Methods METHOD;

    public:
        explicit Bubble(double gamma, double k_rho, double n_int, Methods METHOD = Methods::RUNGE_KUTTA, int i = 0, int j = 0, int k = 0); 
        auto CurveValue(const boost::phoenix::placeholders::arg2_type lambda);
        double LambdaApprox();
        std::string Filename();
        /*  
        double GradVelShock(double gamma_sa, double k_rho);
        double GradVelCon(double gamma_sa, double k_rho);
        */
        double LambdaShockWind(double lambda_c, double lambda_n = 1, double numb_of_years = 1e4, double vel_in = 75000, double r_s = 0);
};
