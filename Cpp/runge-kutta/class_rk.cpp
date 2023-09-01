#include "class_rk.hpp"
#include "const_rk.hpp"
#include <string>
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

double sqr(double x) {
    return x * x;
}

// Observer function
struct StreamingObserver {
    std::ostream &m_out;
    StreamingObserver(std::ostream &out) : m_out(out) {}

    template<class State>
    void operator()(const State &x, double t) const {
        vector_type q = x;
        m_out << t;
        for (size_t i = 0; i < q.size(); ++i) {
            m_out << ' ' << q[i];
        }
        m_out << '\n'; // to get CurveValue(), it should be added here somehow
    }
};

// Value of the curve at lambda
auto Bubble::CurveValue(const boost::phoenix::placeholders::arg2_type lambda) {
    return (gamma__ + 1) / 2 * lambda;
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

        struct System {
            void operator()(const state_type &x, state_type &dxdt, double t) {

                // Derivative of velocity
                dxdt[0] = ((4. * x[0] * gamma_ / (t * (gamma_ + 1.)) - k_rho_ -
                (1. - 1. / eta_) * (x[2] / x[1] * (gamma_ + 1.) / (gamma_ - 1.) * (2. * x[0] / (gamma_ + 1.) - t) * x[0] - 2.))) /
                (x[2] / x[1] * (gamma_ + 1.) / (gamma_ - 1.) * sqr((2. * x[0] / (gamma_ + 1.) - t)) - 2. * gamma_ / (gamma_ + 1.));

                // Derivative of pressure
                dxdt[1] = (-(1 - 1 / eta_) * x[0] - (2 * x[0] / (gamma_ + 1) - t) * dxdt[0]) * (gamma_ + 1) / (gamma_ - 1) * x[2]; // different 

                // initial
                // dxdt[1] = (-(1 - 1 / eta_) * x[0] - (2 * x[0] / (gamma_ + 1) - t) * dxdt[0]) * (gamma_ + 1) / ((gamma_ - 1) * x[2]);

                // Derivative of density
                dxdt[2] = x[2] / gamma_ * ((k_rho_ * (gamma_ - 1) + 2 * (1 - 1 / eta_)) / (2 * x[0] / (gamma_ + 1) - t) + 1 / x[1] * dxdt[1]);
            }
        };
        std::ofstream output;

        output.open(Filename());

        state_type x(3, 1.0); // Size and initial conditions (expecting equal values, but not necessary)
        boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
        
        size_t num_of_steps = integrate_const(stepper, System(), x, x_0, x_c, step,
                output << boost::phoenix::arg_names::arg2 << ' ' << boost::phoenix::arg_names::arg1[0] << ' ' 
                << boost::phoenix::arg_names::arg1[1] << ' ' << boost::phoenix::arg_names::arg1[2] << ' ' << 
                (gamma_ + 1) / 2 * boost::phoenix::arg_names::arg2 << '\n');
        
        // std::clog << num_of_steps << std::endl;
        // std::cout << LambdaApprox() << std::endl;
        output.close();
};

// Approximation using eqn B8a
double Bubble::LambdaApprox() {
    double t = gamma__ * gamma__ * gamma__ + 12 * gamma__ * gamma__ + 8 * gamma__ + 1 - 
    0.5 * (gamma__ + 1) * (3 * gamma__ + 1) * k_rho__ - (gamma__ + 1) * (4 * gamma__ + 1) / eta__;
    double u = 2 * gamma__ * gamma__ * gamma__ + 12 * gamma__ * gamma__ + 
    7 * gamma__ + 1 - 0.5 * (gamma__ + 1) * (3 * gamma__ + 1) * k_rho__ - (gamma__ + 1) * (4 * gamma__ + 1) / eta__;
    return t / u; 
}

std::string Bubble::Filename() {
    std::string g = std::to_string(gamma__);
    boost::trim_right_if(g, boost::is_any_of("0")); // DO NOT CHANGE THE "" to ''
    boost::trim_right_if(g, boost::is_any_of("."));

    std::string k_r = std::to_string(k_rho__);
    boost::trim_right_if(k_r, boost::is_any_of("0"));
    boost::trim_right_if(k_r, boost::is_any_of("."));

    std::string n = std::to_string(n_int__);
    boost::trim_right_if(n, boost::is_any_of("0"));
    boost::trim_right_if(n, boost::is_any_of("."));
    std::string filename = "data/gamma_" + g + "_k_rho_" + k_r + "_n_int_" + n + ".txt";

    return filename;
}