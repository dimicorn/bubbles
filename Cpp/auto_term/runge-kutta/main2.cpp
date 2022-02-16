#include <iostream>
#include "class2.h"
#include <fstream>
#include <vector>

int main(int argc, char **argv) {
    std::vector<double> Gamma, K_rho, N_int;
    /*
    std::cout << "k_rho = 0, n_int = 1" << std::endl;
    // Bubble bubble(5./3., 0, 1, 0, 0, 0); 
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "k_rho = 1, n_int = 1" << std::endl;
    // Bubble bubble1(5./3., 1, 1, 0, 1, 0);
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "k_rho = 2, n_int = 1" << std::endl;
    //Bubble bubble2(5./3., 2, 1, 0, 2, 0);
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "k_rho = 1, n_int = 2" << std::endl;
    // Bubble bubble3(5./3., 1, 2, 0, 1, 1);
    std::cout << "----------------------------------------------" << std::endl;
    std::cout << "k_rho = 0, n_int = 3" << std::endl;
    Bubble bubble4(5./3., 0, 3, 0, 0, 1);
    */

    //Multiple parameters
    for(int i = 0; i < 10; ++i) {
        Gamma.push_back(1.6 + (i + 1) * 0.05);
        K_rho.push_back(0.15 * i);
        N_int.push_back(0.15 * i);
    }

    //OpenMP nested for loop
    #pragma omp parallel for collapse(3)
    for(int i = 0; i < Gamma.size(); ++i) {
        for(int j = 0; j < K_rho.size(); ++j) {
            for (int k = 0; k < N_int.size(); ++k) {
                std::cout << "gamma = " << Gamma[i] << " k_rho = " << K_rho[j] << " n_int =  " << N_int[k] << "\n";
                std::cout << "lambda velocity pressure density" << "\n";
                Bubble bubble(Gamma[i], K_rho[j], N_int[k], i, j, k);
                
            }
        }
    }

    return 0;
}
