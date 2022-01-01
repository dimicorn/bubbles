#include "class.h"
#include "constants.h"
#include <string>
#include <fstream>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

Bubble::Bubble(double gamma, double k_rho, double n_int):
    gamma_(gamma), k_rho_(k_rho), n_int_(n_int) {
        struct stiff_system {
            void operator()(const vector_type &x, vector_type &dxdt, double t) {
                // Derivative of velocity
                dxdt[0] = ((4 * x[0] * gamma_ / (t * (gamma_ + 1)) - k_rho_ - (1 - 1 / eta_) * (x[2] / x[1] * (gamma_ + 1) / (gamma_ - 1) * (2 * x[0] / (gamma_ + 1) - t) * x[0] - 2))) / (x[2] / x[1] * (gamma_ + 1) / (gamma_ - 1) * (2 * x[0] / (gamma_ + 1) - t) * (2 * x[0] / (gamma_ + 1) - t) - 2 * gamma_ / (gamma_ + 1));
                // Derivative of pressure
                dxdt[1] = (-(1 - 1 / eta_) * x[0] - (2 * x[0] / (gamma_ + 1) - t) * dxdt[0]) * (gamma_ + 1) / (gamma_ - 1) * x[2];
                // Derivative of density
                dxdt[2] = x[2] / gamma_ * ((k_rho_ * (gamma_ - 1) + 2 * (1 - 1 / eta_)) / (2 * x[0] / (gamma_ + 1) - t) + 1 / x[1] * dxdt[2]);
            }
        };
        struct stiff_system_jacobi {
            void operator()(const vector_type &x, matrix_type &J, const double &t, vector_type &dfdt) {
                // Jacobian of the system should be changed
                J(0, 0) = -101.0;
                J(0, 1) = -100.0;
                J(1, 0) = 1.0;
                J(1, 1) = 0.0;
                dfdt[0] = 0.0;
                dfdt[1] = 0.0;
            }
        };
        vector_type x(2, 1.0); // Not sure what this is, probably should be changed as well
        size_t num_of_steps = integrate_const(boost::numeric::odeint::make_dense_output<boost::numeric::odeint::rosenbrock4<double>>(1.0e-6, 1.0e-6),
                std::make_pair(stiff_system(), stiff_system_jacobi()),
                x, 0.0, 50.0, 0.01,
                std::cout << boost::phoenix::arg_names::arg2 << " " << boost::phoenix::arg_names::arg1[0] << "\n");
        // std::clog << num_of_steps << std::endl;
    };

// Value of the curve at lambda_c
double Bubble::CurveValue(double lambda_c) {
    return (gamma_ + 1) / 2 * lambda_c;
}

// Approximation using eqn B8a
double Bubble::LambdaApprox() {
    double t = gamma_ * gamma_ * gamma_ + 12 * gamma_ * gamma_ + 8 * gamma_ + 1 - 0.5 * (gamma_ + 1) * (3 * gamma_ + 1) * k_rho_ - (gamma_ + 1) * (4 * gamma_ + 1) / eta_;
    double u = 3 * gamma_ * gamma_ * gamma_ + 12 * gamma_ * gamma_ + 7 * gamma_ + 1 - 0.5 * (gamma_ + 1) * (3 * gamma_ + 1) * k_rho_ - (gamma_ + 1) * (4 * gamma_ + 1) / eta_;
    return t / u; 
}

// Gradient of velocity, r = R_s
double Bubble::GradVel1() {
    double temp = (-(7 * gamma_ + 3) + (gamma_ + 1) * k_rho_ + 3 * (gamma_ + 1) / eta_) / (gamma_ + 1);
    return temp;
}

// Gradient of velocity, r = R_s
double Bubble::GradVel2() {
    double temp = (-2 * (gamma_ + 1) + k_rho_ + 2 / eta_) / (2 * gamma_ / (gamma_ + 1));
    return temp;
}

// Difference between lambda from the ODE and the approximation using eqn B8a
double Bubble::Delta(double lambda_c, double lambda_c_approx) {
    double temp = abs(lambda_c - lambda_c_approx) / lambda_c * 100;
    return temp;
}

// Calculating R_sw
double Bubble::R_sw(double vel, double numb, double k) {
    double f_rho = pow((4 * gamma_ / (gamma_ + 1) * (gamma_ + 1)), 1 / (gamma_ - 1));
    double time = 3600 * 24 * 365 * numb;
    double temp = (gamma_ - 1) / (gamma_ + 1) * vel * f_rho * time * (3 * gamma_ / (3 * (gamma_ - 1) * eta_ + n_int_)) / (k_rho_ * k_rho_ * k_rho_);
    return temp / au;
}

int Bubble::Values() {
    return 0;
}

int Bubble::Norm1() {
    return 0;
}

int Bubble::Norm2() {
    return 0;
}

int Bubble::Q_p() {
     return 0;
}

std::string Bubble::Output() {
    return "Test";
}
