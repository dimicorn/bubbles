#include <iostream>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;


class Bubble{
public:
    double gamma, k_rho, n_int;
    std::vector<double> lambda_c, velocity, density, pressure;
    double q_p(){
        p_r = 2 / (gamma + 1) * pressure[-1];
        p_sw = (2/(gamma + 1)) * pow(((gamma + 1)*(gamma + 1)/(4*gamma)), gamma/(gamma - 1));
        return p_sw/p_r;
    }
private:
    double p_r, p_sw;
};

void rhs( const double x , double &dxdt , const double t )
{
    dxdt = 3.0/(2.0*t*t) + x/(2.0*t);
}

void write_cout( const double &x , const double t )
{
    std::cout << t << '\t' << x << std::endl;
}
typedef runge_kutta_dopri5<double> stepper_type;

int main(){

    std::vector<double> G, K, N_int;
    G = {1.001, 1.01, 1.4, 1.5, 1.6, 5.0/3.0, 1.7, 1.9, 2, 4};
    K = {0, 0.5, 1, 1.5, 2, 2.5, 3};
    N_int = {0, 0.5, 1, 1.5, 2, 2.5, 3};
    for (int i = 0; i < G.size(); ++i){
        for (int j = 0; j < K.size(); ++j){
            for (int q = 0; q < N_int.size(); ++q){
                Bubble bubble;
            }
        }
    }

    double x = 0.0;
    integrate_adaptive(make_controlled(1E-12, 1E-12, stepper_type() ), rhs, x, 1.0, 10.0, 0.1, write_cout);
    return 0;
}
