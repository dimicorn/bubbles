#include "bubble.hpp"

void solve(double gamma, double k_rho, double n_int, Methods m = Methods::RUNGE_KUTTA) {
	Bubble bubble(gamma, k_rho, n_int, m);
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

	std::cout << "lambda_c = " << lambda_c << '\n';
	std::cout << "lambda_sw = " << lambda_sw << '\n';
}
/*
void run(int argc, char** argv) {
	if (argc == 1) {
		// Runge-Kutta
		solve();

	} else if (argc == 2) {
		if (!strcmp(argv[1], "-rk") || !strcmp(argv[1], "--runge-kutta")) {
			// Runge-Kutta
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
*/

int main(int argc, char** argv) {
    int n = 2;
	std::vector<double> Gamma(n), K_rho(n), N_int(n);

    // initial range of gamma, k_rho, n_int
    for (int i = 0; i < n; ++i) {
        Gamma[i] = 2 + 0.1 * i;
        K_rho[i] = 2 + 0.3 * i;
        N_int[i] = 1 + 0.3 * i;
    }

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
	int i_beg, i_end;
    i_beg = rank * n / size;
    if (rank != size - 1) {
        i_end = (rank + 1) * n / size;
    } else {
        i_end = n;
    }
	std::cout << size << '\n';
	std::cout << "I work! rank = " << rank << '\n';
	std::cout << i_beg << ' ' << i_end << '\n';
	
	try {
		if (argc == 1) {
			// Runge-Kutta
			// if (rank == size - 1) {
				for (int i = i_beg; i < i_end; ++i) {
					for (int j = 0; j < n; ++j) {
						for (int k = 0; k < n; ++k) {
							solve(Gamma[i], K_rho[j], N_int[k], Methods::RUNGE_KUTTA);
						}
					}
				}
			/*
			} else {
				for (int i = i_beg; i < i_end; ++i) {
					for (int j = 0; j < n; ++j) {
						for (int k = 0; k < n; ++k) {
							solve(Gamma[i], K_rho[j], N_int[k], Methods::RUNGE_KUTTA);
						}
					}	
				}
			}
			*/

		} else if (argc == 2) {
			if (!strcmp(argv[1], "-rk") || !strcmp(argv[1], "--runge-kutta")) {
				// Runge-Kutta
				// if (rank == size - 1) {
					for (int i = i_beg; i < i_end; ++i) {
						for (int j = 0; j < n; ++j) {
							for (int k = 0; k < n; ++k) {
								solve(Gamma[i], K_rho[j], N_int[k], Methods::RUNGE_KUTTA);
							}
						}	
					}
				/*
				} else {
					for (int i = i_beg; i < i_end; ++i) {
						for (int j = 0; j < n; ++j) {
							for (int k = 0; k < n; ++k) {
								solve(Gamma[i], K_rho[j], N_int[k], Methods::RUNGE_KUTTA);
							}
						}	
					}
				}
				*/
			} else if (!strcmp(argv[1], "-r") || !strcmp(argv[1], "--rosenbrock")) {
				// Rosenbrock
				// if (rank == size - 1) {
					for (int i = i_beg; i < i_end; ++i) {
						for (int j = 0; j < n; ++j) {
							for (int k = 0; k < n; ++k) {
								Bubble bubble(Gamma[i], K_rho[j], N_int[k], Methods::ROSENBROCK);
							}
						}
					}
				/*
				} else {
					for (int i = i_beg; i < i_end; ++i) {
						for (int j = 0; j < n; ++j) {
							for (int k = 0; k < n; ++k) {
								Bubble bubble(Gamma[i], K_rho[j], N_int[k], Methods::ROSENBROCK);
							}
						}
					}
				}
				*/
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
	} catch (std::invalid_argument &e) {
		std::cerr << e.what() << '\n';	
		MPI_Finalize();
		return -1;
	} catch (std::length_error &e) {
		std::cerr << e.what() << '\n';	
		MPI_Finalize();
		return -1;
	}

	MPI_Finalize();

    return 0;
}
