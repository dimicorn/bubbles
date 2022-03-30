#include "class2.h"
#include "constants2.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

typedef std::vector<double> state_type;

double sqr(double x) {
    return x * x;
}

Bubble::Bubble(double gamma, double k_rho, double n_int, int i, int j, int k):
    gamma__(gamma), k_rho__(k_rho), n_int__(n_int), i_(i), j_(j), k_(k) {
        static double gamma_ = gamma__;
        static double k_rho_ = k_rho__;
        static double n_int_ = n_int__;
        static double eta_ = eta__;
        // Not very pretty, but works fine
        if (gamma_ != gamma__) {
            double temp = gamma__ - gamma_;
            gamma_ += temp;
        }
        if (k_rho_ != k_rho__) {
            double temp = k_rho__ - k_rho_;
            k_rho_ += temp;
        }
        if (n_int_ != n_int__) {
            double temp = n_int__ - n_int_;
            n_int_ += temp;
        }
        if (eta_ != eta__) {
            double temp = eta__ - eta_;
            eta_ += temp;
        }

        struct StiffSystem {
            void operator()(const state_type &x, state_type &dxdt, double t) {
                // Derivative of velocity
                dxdt[0] = ((4. * x[0] * gamma_ / (t * (gamma_ + 1.)) - k_rho_ -
                (1. - 1. / eta_) * (x[2] / x[1] * (gamma_ + 1.) / (gamma_ - 1.) * (2. * x[0] / (gamma_ + 1.) - t) * x[0] - 2.))) /
                (x[2] / x[1] * (gamma_ + 1.) / (gamma_ - 1.) * sqr((2. * x[0] / (gamma_ + 1.) - t)) - 2. * gamma_ / (gamma_ + 1.));

                // Derivative of pressure
                dxdt[1] = (-(1 - 1 / eta_) * x[0] - (2 * x[0] / (gamma_ + 1) - t) * dxdt[0]) * (gamma_ + 1) / (gamma_ - 1) * x[2];

                // Derivative of density
                dxdt[2] = x[2] / gamma_ * ((k_rho_ * (gamma_ - 1) + 2 * (1 - 1 / eta_)) / (2 * x[0] / (gamma_ + 1) - t) + 1 / x[1] * dxdt[1]);
            }
        };
        std::ofstream output;
        output.open("data/gamma_" + std::to_string(gamma_) + "_k_rho_" + std::to_string(k_rho_) + "_n_int_" + std::to_string(n_int_) + ".dat");
        output << "gamma =  " << gamma_ << " k_rho = " << k_rho_ << " n_int = " << n_int_ << "\n";
        output << "lambda velocity pressure density\n";
        state_type x(3, 1.0); // Size and initial conditions (expecting equal values, but not necessary)
        boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
        size_t num_of_steps = integrate_const(stepper, StiffSystem(), x,
                1.0, LambdaApprox(), -0.0001,
                output << boost::phoenix::arg_names::arg2 << " " << boost::phoenix::arg_names::arg1[0] << " " << boost::phoenix::arg_names::arg1[1] << " " << boost::phoenix::arg_names::arg1[2] << "\n");
        output.close();
        // std::clog << num_of_steps << std::endl;
        // std::cout << LambdaApprox() << std::endl;
};

// Value of the curve at lambda_c
double Bubble::CurveValue(double lambda_c) {
    return (gamma__ + 1) / 2 * lambda_c;
}

// Approximation using eqn B8a
double Bubble::LambdaApprox() {
    double t = gamma__ * gamma__ * gamma__ + 12 * gamma__ * gamma__ + 8 * gamma__ + 1 - 0.5 * (gamma__ + 1) * (3 * gamma__ + 1) * k_rho__ - (gamma__ + 1) * (4 * gamma__ + 1) / eta__;
    double u = 2 * gamma__ * gamma__ * gamma__ + 12 * gamma__ * gamma__ + 7 * gamma__ + 1 - 0.5 * (gamma__ + 1) * (3 * gamma__ + 1) * k_rho__ - (gamma__ + 1) * (4 * gamma__ + 1) / eta__;
    return t / u; 
}
