#include <iostream>
#include <fstream>
#include "class.h"
#include <vector>

int main(int argc, char **argv) {
    std::vector<float> Gamma, K_rho, N_int;
    // Test parameters
    for (int i = 0; i < 3; ++i) {
        Gamma.push_back(1 + (i + 1) * 0.05);
        K_rho.push_back(0.15 * i);
        N_int.push_back(0.15 * i);
    }
    /*
    // Parallel nested for loop
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < Gamma.size(); ++i) {
        for (int j = 0; j < K_rho.size(); ++j) {
            for (int k = 0; k < N_int.size(); ++k) {
                std::cout << Gamma[i] << ' ' << K_rho[j] << ' ' << N_int[k] << std::endl;
                Bubble bubble(Gamma[i], K_rho[j], N_int[k], i , j , k);
            }
        }
    }
    */
    Bubble bubble(5.0 / 3.0, 0, 3, 0, 0, 0);
    return 0;
}
