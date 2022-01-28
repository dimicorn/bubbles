#include <iostream>
#include <fstream>
#include "class.h"
#include <vector>

int main(int argc, char **argv) {
    std::ofstream output_file;
    std::vector<float> Gamma, K_rho, N_int;
    for (int i = 0; i < 20; ++i) {
        Gamma.push_back(1 + i * 0.05);
        K_rho.push_back(0.15 * i);
        N_int.push_back(0.15 * i);
        std::cout << Gamma[i] << ' ' << K_rho[i] << ' ' << N_int[i];
    }
    output_file.open("output.txt");
    for (int i = 0; i < Gamma.size(); ++i) {
        for (int j = 0; j < K_rho.size(); ++j) {
            for (int k = 0; k < N_int.size(); ++k) {
                Bubble bubble(Gamma[i], K_rho[j], N_int[k]);
                output_file << bubble.Output() << std::endl;
            }
        }
    }
    return 0;
}
