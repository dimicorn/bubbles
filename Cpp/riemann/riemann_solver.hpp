#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

class RiemannSolver {
private:
    double gamma, g1, g2, g3, g4, g5, g6, g7, g8;
    double dl, ul, pl, cl, dr, ur, pr, cr;

    std::string output_file = "exact.out";

    void starpu(double &p, double &u, const double &mpa);
    void guessp(double &pm);
    void prefun(double &f, double &fd, double &p, double &dk, double &pk, double &ck);
    void sample(double &pm, double &um, double &s, double &d, double &u, double &p);

public:
    RiemannSolver(const double &domlen, const double &diaph, const int &cells, 
    const double &gamma, const double &timeout, const double &dl, const double &ul, 
    const double &pl, const double &dr, const double &ur, const double &pr, 
    const double &mpa);
};
