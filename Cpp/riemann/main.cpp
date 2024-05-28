#include "riemann_solver.hpp"

int main(int argc, char** argv) {
    if (argc == 1) {
        std::cout << "No input file!\n";
        std::cout << "Example: g++ riemann.cpp input_file\n";
        return 0;
    }
    int cells;
    double domlen, diaph, gamma, timeout;
    double dl, ul, pl, dr, ur, pr, mpa;

    // TESTING 
    /*
    std::string test_file = "test.ini";
    std::ifstream test;
    test.open(test_file);
    for (int i = 0; i < 2; ++i) {
        test >> domlen >> diaph >> cells >> gamma >> timeout;
        test >> dl >> ul >> pl >> dr >> ur >> pr >> mpa;
        RiemannSolver rs(domlen, diaph, cells, gamma, timeout, dl, ul, pl, dr, ur, pr, mpa);
    }
    test.close();
    */

    std::string input_file = argv[1];
    std::ifstream input;
    input.open(input_file);
    input >> domlen >> diaph >> cells >> gamma >> timeout;
    input >> dl >> ul >> pl >> dr >> ur >> pr >> mpa;
    input.close();

    RiemannSolver rs(domlen, diaph, cells, gamma, timeout, dl, ul, pl, dr, ur, pr, mpa);

    return 0;
}
