#include <iostream>
#include "class_rk.hpp"
#include <fstream>
#include <vector>
#include "mpi.h"


void vectorize(std::vector<double> v, double beg, double end, double step) {
    for (int i = beg; i <= end; i += step) {
        v.push_back(i);
    }
}

int main(int argc, char** argv) {
    std::vector<double> Gamma, K_rho, N_int;
    
    // initial range of gamma, k_rho, n_int
    int n = 12;
    for (int i = 0; i <= n; ++i) {
        Gamma.push_back(1.2 + 0.1 * i);
        K_rho.push_back(0.3 * i);
        N_int.push_back(0.3 * i);
    }

    /*
    std::ifstream input;
    input.open("input.dat");
    if (!input) {
        std::cout << "Input file couldn't be opened!\n";
    }

    double gamma_beg, gamma_end, gamma_step;
    input >> gamma_beg >> gamma_end >> gamma_step;
    vectorize(Gamma, gamma_beg, gamma_end, gamma_step);
    int n = Gamma.size();

    double k_rho_beg, k_rho_end, k_rho_step;
    input >> k_rho_beg >> k_rho_end >> k_rho_step;
    vectorize(K_rho, k_rho_beg, k_rho_end, k_rho_step);

    double n_int_beg, n_int_end, n_int_step;
    input >> n_int_beg >> n_int_end >> n_int_step;
    vectorize(N_int, n_int_beg, n_int_end, n_int_step);
     
    input.close();
    */

    int rank, size;
    
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int i_beg, i_end;
    i_beg = rank * n / size;
    if (rank != size - 1) {
        i_end = (rank + 1) * n / size;
    } else {
        i_end = n;
    }
    std::cout << "I work!" << rank << '\n';
    if (rank == size - 1) {
        for (int i = i_beg; i <= i_end; ++i) {
            for (int j = 0; j < K_rho.size(); ++j) {
                for (int k = 0; k < N_int.size(); ++k) {
                    Bubble bubble(Gamma[i], K_rho[j], N_int[k], i, j, k);
                }
            }
        }
    } else {
        for (int i = i_beg; i < i_end; ++i) {
            for (int j = 0; j < K_rho.size(); ++j) {
                for (int k = 0; k < N_int.size(); ++k) {
                    Bubble bubble(Gamma[i], K_rho[j], N_int[k], i, j, k);
                }
            }
        }
    }

    MPI_Finalize();

    return 0;
}
