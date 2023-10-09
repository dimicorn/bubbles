#include "bubble.hpp"

const std::string unknown_arg = "unknown argument: ";
const std::string too_many_args = "too many arguments";

void solve() {
	Bubble bubble(2, 2, 1, Methods::RUNGE_KUTTA);
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
}

void run(int argc, char** argv) {
	if (argc == 1) {
		// Runge-Kutta
		solve();

	} else if (argc == 2) {
		if (!strcmp(argv[1], "-rk") || !strcmp(argv[1], "--runge-kutta")) {
			// Runge-Kutta
			// METHOD = Methods::RUNGE_KUTTA;
			solve();
		} else if (!strcmp(argv[1], "-r") || !strcmp(argv[1], "--rosenbrock")) {
			// Rosenbrock
			Bubble bubble(2, 2, 1, Methods::ROSENBROCK);

		} else if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
			// Help
		} else {
			// Throw invalid argument exception
			throw std::invalid_argument(unknown_arg + "'" + argv[1] + "'");
		}
	} else {
		// Throw length error exception
		throw std::length_error(too_many_args);
	}

}


int main(int argc, char** argv) {
	/*
	int size, rank;
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
	*/
	
	try {
		run(argc, argv);
	} catch (std::invalid_argument &e) {
		std::cerr << e.what() << '\n';	
		return -1;
	} catch (std::length_error &e) {
		std::cerr << e.what() << '\n';	
		return -1;
	}

    return 0;
}
