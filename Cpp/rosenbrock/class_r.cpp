#include "class_r.hpp"
#include "const_r.hpp"
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <utility>

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/algorithm/string.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

// Observer function
struct StreamingObserver {
    std::ostream& m_out;
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
            void operator()(const vector_type &x, vector_type &dxdt, double t) {

                // Derivative of velocity
                dxdt[0] = ((4 * x[0] * gamma_ / (t * (gamma_ + 1)) - k_rho_ - (1 - 1 / eta_) * (x[2] / x[1] * (gamma_ + 1) / 
                (gamma_ - 1) * (2 * x[0] / (gamma_ + 1) - t) * x[0] - 2))) / (x[2] / x[1] * (gamma_ + 1) / (gamma_ - 1) * 
                (2 * x[0] / (gamma_ + 1) - t) * (2 * x[0] / (gamma_ + 1) - t) - 2 * gamma_ / (gamma_ + 1));

                // Derivative of pressure
                dxdt[1] = (-(1 - 1 / eta_) * x[0] - (2 * x[0] / (gamma_ + 1) - t) * dxdt[0]) * (gamma_ + 1) / (gamma_ - 1) * x[2]; // different

                // Derivative of density
                dxdt[2] = x[2] / gamma_ * ((k_rho_ * (gamma_ - 1) + 2 * (1 - 1 / eta_)) / (2 * x[0] / (gamma_ + 1) - t) + 1 / x[1] * dxdt[1]);
            }
        };
        struct StiffSystemJacobi {
            void operator()(const vector_type &x, matrix_type &J, const double &t, vector_type &dfdt) {
                double alpha = ((4 * x[0] * gamma_ / (t * (gamma_ + 1)) - k_rho_ - (1 - 1 / eta_) * (x[2] / x[1] * (gamma_ + 1) / (gamma_ - 1) * (2 * x[0] / (gamma_ + 1) - t) * x[0] - 2)));
                double beta = (x[2] / x[1] * (gamma_ + 1) / (gamma_ - 1) * (2 * x[0] / (gamma_ + 1) - t) * (2 * x[0] / (gamma_ + 1) - t) - 2 * gamma_ / (gamma_ + 1));

                // Derivatives of velocity, pressure and density
                double dvdl = alpha / beta;
                double dpdl = (-(1 - 1 / eta_) * x[0] - (2 * x[0] / (gamma_ + 1) - t) * dvdl) * (gamma_ + 1) / (gamma_ - 1) * x[2];
                double drhodl = x[2] / gamma_ * ((k_rho_ * (gamma_ - 1) + 2 * (1 - 1 / eta_)) / (2 * x[0] / (gamma_ + 1) - t) + 1 / x[1] * dpdl);

                // Jacobian
                J(0, 0) = ((4 * gamma_/ (t * (gamma_ + 1)) - (1 - 1 / eta_) * x[2] / x[1] * (gamma_ + 1) / (gamma_ - 1) * (4 * x[0] / (gamma_ + 1) - t)) * beta - alpha * (x[2] / x[1] * (gamma_ + 1) / (gamma_ - 1) * (2 * x[0] / (gamma_ + 1) - t) * 4 / (gamma_ + 1))) / (beta * beta);
                J(0, 1) = ((1 - 1 / eta_) * (-x[2] / (x[1] * x[1]) * (gamma_ + 1) / (gamma_ - 1) * (2 * x[0] / (gamma_ + 1) - t) * x[0]) * beta - alpha * ((4 * x[0] * gamma_ / (t * (gamma_ + 1)) - k_rho_ - (1 - 1 / eta_) * (x[2] / x[1] * (gamma_ + 1) / (gamma_ - 1) * (2 * x[0] / (gamma_ + 1) - t) * x[0] - 2)))) / (beta * beta);
                J(0, 2) = ((1 - 1 / eta_) * 1 / x[1] * (gamma_ + 1) / (gamma_ - 1) * (2 * x[0] / (gamma_ + 1) - t) * x[0] * beta - alpha * ((4 * x[0] * gamma_ / (t * (gamma_ + 1)) - k_rho_ - (1 - 1 / eta_) * (x[2] / x[1] * (gamma_ + 1) / (gamma_ - 1) * (2 * x[0] / (gamma_ + 1) - t) * x[0] - 2)))) / (beta * beta);
                J(1, 0) = (gamma_ + 1) / (gamma_ - 1) * x[2] * (-1 + 1 / eta_ - (2 / (gamma_ + 1) * dvdl + (2 * x[0] / (gamma_ + 1) - t) * J(0, 0)));
                J(1, 1) = (gamma_ + 1) / (gamma_ - 1) * -x[2] * (2 * x[0] / (gamma_ + 1) - t) * J(0, 1);
                J(1, 2) = (gamma_ + 1) / (gamma_ - 1) * (-x[0] * (1 - 1 / eta_) - (2 * x[0] / (gamma_ + 1) - t) * dvdl) + (gamma_ + 1) / (gamma_ - 1) * -x[2] * (2 * x[0] / (gamma_ + 1) - t) * J(0, 2);
                J(2, 0) = x[2] / gamma_ * (k_rho_ * (gamma_ - 1) + 2 * (1 - 1 / eta_) * (-2) / ((gamma_ + 1) * (2 * x[0] / (gamma_ + 1) - t) * (2 * x[0] / (gamma_ + 1) - t)) + 1 / x[1] * J(1, 0));
                J(2, 1) = x[2] / gamma_ * (-1 / (x[1] * x[1]) * dpdl + 1 / x[1] * J(1, 1));
                J(2, 2) = 1 / gamma_ * ((k_rho_ * (gamma_ - 1) + 2 * (1 - 1 / eta_)) / (2 * x[0] / (gamma_ + 1) - t) + 1 / x[1] * dpdl) + x[2] / gamma_ * 1 / x[1] * J(1, 2);

                // Derivatives by lambda
                dfdt[0] = ((4 * gamma_ / (gamma_ + 1) * (dvdl / t - x[0] / (t * t)) - (1 - 1 / eta_) * ((gamma_ + 1) / (gamma_ - 1) * (((dpdl * x[0] + dvdl * x[1]) * x[2] - x[0] * x[1] * drhodl) / (x[2] * x[2]) * (2 * x[0] / (gamma_ + 1) - t) + x[0] * x[1] / x[2] * (2 / (gamma_ + 1) * dvdl - 1)))) * beta - alpha * ((gamma_ + 1) / (gamma_ - 1) * ((drhodl * x[1] + dpdl * x[2]) * ((2 * x[0] / (gamma_ + 1) - t) * (2 * x[0] / (gamma_ + 1) - t)) + 2 * x[2] / x[1] * (2 * x[0] / (gamma_ + 1) - t) * (2 / (gamma_ + 1) * dvdl - 1)))) / (beta * beta);
                dfdt[1] = (gamma_ + 1) / (gamma_ - 1) * ((-(1 - 1 / eta_) * dvdl - ((2 / (gamma_ + 1) * dvdl - 1) * dvdl + (2 * x[0] / (gamma_ + 1) - t) * dfdt[0]))+ (-(1 - 1 / eta_) * x[0] - (2 * x[0] / (gamma_ + 1) - t) * dvdl) * drhodl);
                dfdt[2] = 1 / gamma_ * (drhodl * ((k_rho_ * (gamma_ - 1) + 2 * (1 - 1 / eta_)) / (2 * x[0] / (gamma_ + 1) - t) + dpdl / x[1]) + x[2] * ((-k_rho_ * (gamma_ - 1) - 2 * (1 - 1 / eta_)) * (2 / (gamma_ + 1) * dvdl - 1) / ((2 * x[0] / (gamma_ + 1) - t) * (2 * x[0] / (gamma_ + 1) - t)) + (-dpdl / (x[1] * x[1]) * dpdl + dfdt[1] / x[1])));
            }
        };

        std::ofstream output;
        // std::cout << gamma_ << " " << k_rho_ << " " << n_int_ << std::endl;
        
        std::string g = std::to_string(gamma__);
        boost::trim_right_if(g, boost::is_any_of("0")); // DO NOT CHANGE THE "" -> ''
        boost::trim_right_if(g, boost::is_any_of("."));

        std::string k_r = std::to_string(k_rho__);
        boost::trim_right_if(k_r, boost::is_any_of("0"));
        boost::trim_right_if(k_r, boost::is_any_of("."));

        std::string n = std::to_string(n_int__);
        boost::trim_right_if(n, boost::is_any_of("0"));
        boost::trim_right_if(n, boost::is_any_of("."));
        std::string filename = "data/gamma_" + g + "_k_rho_" + k_r + "_n_int_" + n + ".txt";
        output.open(filename);

        vector_type x(3, 1.0); // Size and initial conditions (expecting equal values)

        size_t num_of_steps = integrate_const(boost::numeric::odeint::make_dense_output<boost::numeric::odeint::rosenbrock4<double>>(eps, eps),
                std::make_pair(StiffSystem(), StiffSystemJacobi()),
                x, x_0, LambdaApprox(), step, 
                output << boost::phoenix::arg_names::arg2 << ' ' << boost::phoenix::arg_names::arg1[0] << ' ' 
                << boost::phoenix::arg_names::arg1[1] << ' ' << boost::phoenix::arg_names::arg1[2] << ' '
                << (gamma_ + 1) / 2 * boost::phoenix::arg_names::arg2 << '\n');
                // StreamingObserver(std::cout)); // to make output to file change std::cout to output
        
        // output << boost::phoenix::arg_names::arg2 << ' ' << boost::phoenix::arg_names::arg1[0] << ' ' << boost::phoenix::arg_names::arg1[1] << ' ' << boost::phoenix::arg_names::arg1[2] << '\n');
        // std::clog << num_of_steps << std::endl;
        // std::cout << LambdaApprox() << std::endl;

        output.close();
    };

// Value of the curve at lambda_c
double Bubble::CurveValue(double lambda) {
    return (gamma__ + 1) / 2 * lambda;
}

// Approximation using eqn B8a
double Bubble::LambdaApprox() {
    double t = gamma__ * gamma__ * gamma__ + 12 * gamma__ * gamma__ + 8 * gamma__ + 1 - 0.5 * (gamma__ + 1) * (3 * gamma__ + 1) * k_rho__ - (gamma__ + 1) * (4 * gamma__ + 1) / eta__;
    double u = 2 * gamma__ * gamma__ * gamma__ + 12 * gamma__ * gamma__ + 7 * gamma__ + 1 - 0.5 * (gamma__ + 1) * (3 * gamma__ + 1) * k_rho__ - (gamma__ + 1) * (4 * gamma__ + 1) / eta__;
    return t / u; 
}

/*
// Gradient of velocity, r = R_s
double Bubble::GradVel1() {
    double temp = (-(7 * gamma__ + 3) + (gamma__ + 1) * k_rho__ + 3 * (gamma__ + 1) / eta__) / (gamma__ + 1);
    return temp;
}

// Gradient of velocity, r = R_s
double Bubble::GradVel2() {
    double temp = (-2 * (gamma__ + 1) + k_rho__ + 2 / eta__) / (2 * gamma__ / (gamma__ + 1));
    return temp;
}

// Difference between lambda from the ODE and the approximation using eqn B8a
double Bubble::Delta(double lambda_c, double lambda_c_approx) {
    double temp = abs(lambda_c - lambda_c_approx) / lambda_c * 100;
    return temp;
}

// Calculating R_sw
double Bubble::R_sw(double vel, double numb, double k) {
    double f_rho = pow((4 * gamma__ / (gamma__ + 1) * (gamma__ + 1)), 1 / (gamma__ - 1));
    double time = 3600 * 24 * 365 * numb;
    double temp = (gamma__ - 1) / (gamma__ + 1) * vel * f_rho * time * (3 * gamma__ / (3 * (gamma__ - 1) * eta__ + n_int__)) / (k_rho__ * k_rho__ * k_rho__);
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

void Bubble::Observer() {
    std::cout << "hi\n" << std::endl;
}
*/