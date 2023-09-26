#include "class_rk.hpp"

int main(int argc, char** argv) {
    Bubble bubble(1.4, 2, 1);
	double lambda_n = 0.9;

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

	double lambda_sw = bubble.LambdaShockWind(lambda_c, lambda_n);

	printf("lambda_c = %lf\n", lambda_c);
	printf("lambda_sw = %lf\n", lambda_sw);

    return 0;
}
