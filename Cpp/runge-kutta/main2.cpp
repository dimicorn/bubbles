#include <iostream>
#include "class2.h"
#include <fstream>
#include <vector>

int main(int argc, char **argv) {
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
    std::cout << "k_rho = 0, n_int = 3" << std::endl;
    Bubble bubble4(5./3., 0, 3, 0, 0, 1);
    return 0;
}
