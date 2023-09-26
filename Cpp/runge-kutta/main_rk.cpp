#include "class_rk.hpp"

int main(int argc, char** argv) {
	/*
    std::cout << "k_rho = 0, n_int = 1" << std::endl;
    Bubble bubble(5./3., 0, 1, 0, 0, 0); 
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "k_rho = 1, n_int = 1" << std::endl;
    Bubble bubble1(5./3., 1, 1, 0, 1, 0);
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "k_rho = 2, n_int = 1" << std::endl;
    Bubble bubble2(5./3., 2, 1, 0, 2, 0);
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "k_rho = 1, n_int = 2" << std::endl;
    Bubble bubble3(5./3., 1, 2, 0, 1, 1);
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "gamma = 2, k_rho = 0, n_int = 3" << std::endl;
	*/

    Bubble bubble(2.5, 2, 1);

	std::string output_filename = bubble.Filename();

	std::ifstream output;
    output.open(output_filename);
    double lambda, vel, pres, dens, curve;
    double lambda_c, vel_c = 1, curve_c = 0;
    while (output >> lambda >> vel >> pres >> dens >> curve) {
        if (abs(vel - curve) < abs(vel_c - curve_c)) {
            lambda_c = lambda;
            vel_c = vel;
            curve_c = curve;
        }
    }
    output.close();

	std::cout << bubble.LambdaShockWind(lambda_c) << '\n';

    return 0;
}
