#include<iostream>

/* defining ordinary differential equation to be solved */
/* In this example we are solving dy/dx = x + y */
#define f(x,y) x+y

using namespace std;

int main()
{
    float x0, y0, xn, h, yn, slope;
    int i, n;

    cout<<"Enter Initial Condition"<< endl;
    cout<<"x0 = ";
    cin>> x0;
    cout<<"y0 = ";
    cin >> y0;
    cout<<"Enter calculation point xn = ";
    cin>>xn;
    cout<<"Enter number of steps: ";
    cin>> n;

    /* Calculating step size (h) */
    h = (xn-x0)/n;

    /* Euler's Method */
    cout<<"\nx0\ty0\tslope\tyn\n";
    cout<<"------------------------------\n";

    for(i=0; i < n; i++)
    {
        slope = f(x0, y0);
        yn = y0 + h * slope;
        cout<< x0<<"\t"<< y0<<"\t"<< slope<<"\t"<< yn<< endl;
        y0 = yn;
        x0 = x0+h;
    }

    /* Displaying result */
    cout<<"\nValue of y at x = "<< xn<< " is " << yn;

    return 0;
}



/*
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

//typedef runge_kutta_dopri5<double> stepper_type;

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
*/
