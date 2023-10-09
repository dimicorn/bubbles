#include "bubble.hpp"

double sqr(const double &x) {
    return x * x;
}

// Observer function, not used for output
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
    return (gamma + 1) / 2 * lambda;
}

Bubble::Bubble(double gamma, double k_rho, double n_int, Methods METHOD, int i, int j, int k):
gamma(gamma), k_rho(k_rho), n_int(n_int), i(i), j(j), k(k), eta((2 + n_int) / (5 - k_rho)), METHOD(METHOD) {
    static double _gamma = 0, _k_rho = 0, _n_int = 0, _eta = 0;

    // A bit prettier
	_gamma += gamma;
	_k_rho += k_rho;
	_n_int += n_int;
	_eta += eta;

    struct System {
        void operator()(const state_type &x, state_type &dxdt, double t) {

            // Derivative of velocity
            dxdt[0] = ((4. * x[0] * _gamma / (t * (_gamma + 1.)) - _k_rho -
            (1. - 1. / _eta) * (x[2] / x[1] * (_gamma + 1.) / (_gamma - 1.) * (2. * x[0] / (_gamma + 1.) - t) * x[0] - 2.))) /
            (x[2] / x[1] * (_gamma + 1.) / (_gamma - 1.) * sqr((2. * x[0] / (_gamma + 1.) - t)) - 2. * _gamma / (_gamma + 1.));

            // Derivative of pressure
            dxdt[1] = (-(1 - 1 / _eta) * x[0] - (2 * x[0] / (_gamma + 1) - t) * dxdt[0]) * (_gamma + 1) / (_gamma - 1) * x[2]; // different 

            // initial
            // dxdt[1] = (-(1 - 1 / eta) * x[0] - (2 * x[0] / (gamma + 1) - t) * dxdt[0]) * (gamma + 1) / ((gamma - 1) * x[2]);

            // Derivative of density
            dxdt[2] = x[2] / _gamma * ((_k_rho * (_gamma - 1) + 2 * (1 - 1 / _eta)) / (2 * x[0] / (_gamma + 1) - t) + 1 / x[1] * dxdt[1]);
    	}
    };

	struct StiffSystem {
		void operator()(const vector_type &x, vector_type &dxdt, double t) {

            // Derivative of velocity
            dxdt[0] = ((4 * x[0] * _gamma / (t * (_gamma + 1)) - _k_rho - (1 - 1 / _eta) * (x[2] / x[1] * (_gamma + 1) /
			(_gamma - 1) * (2 * x[0] / (_gamma + 1) - t) * x[0] - 2))) / (x[2] / x[1] * (_gamma + 1) / (_gamma - 1) *
            (2 * x[0] / (_gamma + 1) - t) * (2 * x[0] / (_gamma + 1) - t) - 2 * _gamma / (_gamma + 1));

            // Derivative of pressure
            dxdt[1] = (-(1 - 1 / _eta) * x[0] - (2 * x[0] / (_gamma + 1) - t) * dxdt[0]) * (_gamma + 1) / (_gamma - 1) * x[2]; // different

            // Derivative of density
            dxdt[2] = x[2] / _gamma * ((_k_rho * (_gamma - 1) + 2 * (1 - 1 / _eta)) / (2 * x[0] / (_gamma + 1) - t) + 1 / x[1] * dxdt[1]);
        }
    }; 

	struct StiffSystemJacobi {
        void operator()(const vector_type &x, matrix_type &J, const double &t, vector_type &dfdt) {
            double alpha = ((4 * x[0] * _gamma / (t * (_gamma + 1)) - _k_rho - (1 - 1 / _eta) * (x[2] / x[1] * (_gamma + 1) / (_gamma - 1) * (2 * x[0] / (_gamma + 1) - t) * x[0] - 2)));
            double beta = (x[2] / x[1] * (_gamma + 1) / (_gamma - 1) * (2 * x[0] / (_gamma + 1) - t) * (2 * x[0] / (_gamma + 1) - t) - 2 * _gamma / (_gamma + 1));

            // Derivatives of velocity, pressure and density
            double dvdl = alpha / beta;
            double dpdl = (-(1 - 1 / _eta) * x[0] - (2 * x[0] / (_gamma + 1) - t) * dvdl) * (_gamma + 1) / (_gamma - 1) * x[2];
            double drhodl = x[2] / _gamma * ((_k_rho * (_gamma - 1) + 2 * (1 - 1 / _eta)) / (2 * x[0] / (_gamma + 1) - t) + 1 / x[1] * dpdl);

            // Jacobian
            J(0, 0) = ((4 * _gamma/ (t * (_gamma + 1)) - (1 - 1 / _eta) * x[2] / x[1] * (_gamma + 1) / (_gamma - 1) * (4 * x[0] / (_gamma + 1) - t)) * beta -
			alpha * (x[2] / x[1] * (_gamma + 1) / (_gamma - 1) * (2 * x[0] / (_gamma + 1) - t) * 4 / (_gamma + 1))) / (beta * beta);
            J(0, 1) = ((1 - 1 / _eta) * (-x[2] / (x[1] * x[1]) * (_gamma + 1) / (_gamma - 1) * (2 * x[0] / (_gamma + 1) - t) * x[0]) * beta -
			alpha * ((4 * x[0] * _gamma / (t * (_gamma + 1)) - _k_rho - (1 - 1 / _eta) * (x[2] / x[1] * (_gamma + 1) / (_gamma - 1) * (2 * x[0] / (_gamma + 1) - t) * x[0] - 2)))) / (beta * beta);
            J(0, 2) = ((1 - 1 / _eta) * 1 / x[1] * (_gamma + 1) / (_gamma - 1) * (2 * x[0] / (_gamma + 1) - t) * x[0] * beta -
			alpha * ((4 * x[0] * _gamma / (t * (_gamma + 1)) - _k_rho - (1 - 1 / _eta) * (x[2] / x[1] * (_gamma + 1) / (_gamma - 1) * (2 * x[0] / (_gamma + 1) - t) * x[0] - 2)))) / (beta * beta);
            J(1, 0) = (_gamma + 1) / (_gamma - 1) * x[2] * (-1 + 1 / _eta - (2 / (_gamma + 1) * dvdl + (2 * x[0] / (_gamma + 1) - t) * J(0, 0)));
            J(1, 1) = (_gamma + 1) / (_gamma - 1) * -x[2] * (2 * x[0] / (_gamma + 1) - t) * J(0, 1);
            J(1, 2) = (_gamma + 1) / (_gamma - 1) * (-x[0] * (1 - 1 / _eta) - (2 * x[0] / (_gamma + 1) - t) * dvdl) + (_gamma + 1) / (_gamma - 1) * -x[2] * (2 * x[0] / (_gamma + 1) - t) * J(0, 2);
            J(2, 0) = x[2] / _gamma * (_k_rho * (_gamma - 1) + 2 * (1 - 1 / _eta) * (-2) / ((_gamma + 1) * (2 * x[0] / (_gamma + 1) - t) * (2 * x[0] / (_gamma + 1) - t)) + 1 / x[1] * J(1, 0));
            J(2, 1) = x[2] / _gamma * (-1 / (x[1] * x[1]) * dpdl + 1 / x[1] * J(1, 1));
            J(2, 2) = 1 / _gamma * ((_k_rho * (_gamma - 1) + 2 * (1 - 1 / _eta)) / (2 * x[0] / (_gamma + 1) - t) + 1 / x[1] * dpdl) + x[2] / _gamma * 1 / x[1] * J(1, 2);

            // Derivatives by lambda
            dfdt[0] = ((4 * _gamma / (_gamma + 1) * (dvdl / t - x[0] / (t * t)) - (1 - 1 / _eta) * ((_gamma + 1) / (_gamma - 1) * 
			(((dpdl * x[0] + dvdl * x[1]) * x[2] - x[0] * x[1] * drhodl) / (x[2] * x[2]) * (2 * x[0] / (_gamma + 1) - t) + x[0] * x[1] / x[2] * (2 / (_gamma + 1) * dvdl - 1)))) * 
			beta - alpha * ((_gamma + 1) / (_gamma - 1) * ((drhodl * x[1] + dpdl * x[2]) * ((2 * x[0] / (_gamma + 1) - t) * (2 * x[0] / (_gamma + 1) - t)) + 
			2 * x[2] / x[1] * (2 * x[0] / (_gamma + 1) - t) * (2 / (_gamma + 1) * dvdl - 1)))) / (beta * beta);

            dfdt[1] = (_gamma + 1) / (_gamma - 1) * ((-(1 - 1 / _eta) * dvdl - ((2 / (_gamma + 1) * dvdl - 1) * dvdl + (2 * x[0] / (_gamma + 1) - t) * dfdt[0])) +
			(-(1 - 1 / _eta) * x[0] - (2 * x[0] / (_gamma + 1) - t) * dvdl) * drhodl);

            dfdt[2] = 1 / _gamma * (drhodl * ((_k_rho * (_gamma - 1) + 2 * (1 - 1 / _eta)) / (2 * x[0] / (_gamma + 1) - t) + dpdl / x[1]) + 
			x[2] * ((-_k_rho * (_gamma - 1) - 2 * (1 - 1 / _eta)) * (2 / (_gamma + 1) * dvdl - 1) / ((2 * x[0] / (_gamma + 1) - t) * (2 * x[0] / (_gamma + 1) - t)) +
			(-dpdl / (x[1] * x[1]) * dpdl + dfdt[1] / x[1])));
            }
        };

    std::ofstream output;
    output.open(Filename());

	if (METHOD == Methods::RUNGE_KUTTA) {
    	state_type x(3, 1.0); // Size and initial conditions (expecting equal values, but not necessary)
    	boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
        
    	size_t num_of_steps = integrate_const(stepper, System(), x, x_0, x_c, step,
				output << boost::phoenix::arg_names::arg2 << ' ' << boost::phoenix::arg_names::arg1[0] << ' ' 
            	<< boost::phoenix::arg_names::arg1[1] << ' ' << boost::phoenix::arg_names::arg1[2] << ' ' << 
            	(gamma + 1) / 2 * boost::phoenix::arg_names::arg2 << '\n');
        
    	// std::clog << num_of_steps << std::endl;
    	// std::cout << LambdaApprox() << std::endl;
	} else if (METHOD == Methods::ROSENBROCK) {
		vector_type x(3, 1.0); // Size and initial conditions (expecting equal values)

        // const boost::phoenix::placeholders::arg1_type Y;
        // const boost::phoenix::placeholders::arg2_type X;

        size_t num_of_steps = integrate_const(boost::numeric::odeint::make_dense_output<boost::numeric::odeint::rosenbrock4<double>>(eps, eps),
                std::make_pair(StiffSystem(), StiffSystemJacobi()),
                x, x_0, LambdaApprox(), step,
                output << boost::phoenix::arg_names::arg2 << ' ' << boost::phoenix::arg_names::arg1[0] << ' '
                << boost::phoenix::arg_names::arg1[1] << ' ' << boost::phoenix::arg_names::arg1[2] << ' '
                << CurveValue(boost::phoenix::arg_names::arg2) << '\n');
                // StreamingObserver(std::cout)); // to make output to file change std::cout to output

        // std::cout << LambdaApprox() << std::endl;
	}
	output.close();
};

// Approximation using eqn B8a
double Bubble::LambdaApprox() {
    double t = sqr(gamma) * gamma + 12 * sqr(gamma) + 8 * gamma + 1 - 
    0.5 * (gamma + 1) * (3 * gamma + 1) * k_rho - (gamma + 1) * (4 * gamma + 1) / eta;
    double u = 2 * sqr(gamma) * gamma + 12 * sqr(gamma) + 
    7 * gamma + 1 - 0.5 * (gamma + 1) * (3 * gamma + 1) * k_rho - (gamma + 1) * (4 * gamma + 1) / eta;
    return t / u; 
}

std::string Bubble::Filename() {
    std::string g = std::to_string(gamma);
    boost::trim_right_if(g, boost::is_any_of("0")); // DO NOT CHANGE THE "" to ''
    boost::trim_right_if(g, boost::is_any_of("."));

    std::string k_r = std::to_string(k_rho);
    boost::trim_right_if(k_r, boost::is_any_of("0"));
    boost::trim_right_if(k_r, boost::is_any_of("."));

    std::string n = std::to_string(n_int);
    boost::trim_right_if(n, boost::is_any_of("0"));
    boost::trim_right_if(n, boost::is_any_of("."));
    std::string filename = "data/gamma_" + g + "_k_rho_" + k_r + "_n_int_" + n + ".txt";

    return filename;
}

// Appendix B2
// lambda_sw = R_s / R_sw (shocked wind)
double Bubble::LambdaShockWind(double lambda_c, double lambda_n, double numb_of_years, double vel_in, double r_s) {
	double gamma_sw = gamma; // For now it is true

    double f_rho = pow((4 * gamma_sw / sqr(gamma_sw + 1)), 1 / (gamma_sw - 1));
    double time = numb_of_years * year_in_sec;
	// 1.165 *
	double lambda_sw2 = lambda_n * (gamma_sw - 1) / (gamma_sw + 1) * f_rho / (sqr(lambda_c) * lambda_c * (3 * (gamma_sw -  1) * eta + n_int) / (3 * gamma_sw));

	return sqrt(lambda_sw2);

	/*
    double numerator = r_s * sqr(r_s) * lambda_c * sqr(lambda_c) * (3 * (gamma_sw -  1) * eta + n_int) / (3 * gamma_sw);
    double r_sw = sqrt(numerator / ((gamma_sw - 1) / (gamma_sw + 1) * vel_in * f_rho * time));
    return r_sw / au; // in astronomical units
	*/
}

// n_int = 1, independent of wind
// t ~ 10^4(3) - 10^5 years
// v_in ~ 50-100 km/s
// R_s ~ v_in * time
// 0.1 - 2 => 0.5 - 1
// R_s / R_sw > 1 => R_c / R_sw > 1

/*
// Gradient of velocity, r = R_s (shocked)
double Bubble::GradVelShock(double gammasa, double k_rho) {
    return (-(7 * gammasa + 3) + (gammasa + 1) * k_rho + 3 * (gammasa + 1) / eta) / (gammasa + 1);
}

// Gradient of velocity, r = R_c (at discontinuity)
double Bubble::GradVelCon(double gammasa, double k_rho) {
    return (-2 * (gammasa + 1) + k_rho + 2 / eta) / (2 * gammasa / (gammasa + 1));
}
*/
