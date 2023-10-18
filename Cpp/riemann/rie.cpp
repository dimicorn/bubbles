#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

class RiemannSolver {
private:
    int cells;
    double gamma, g1, g2, g3, g4, g5, g6, g7, g8;
    double dl, ul, pl, cl, dr, ur, pr, cr;

    double diaph, domlen, ds, dx, pm, mpa, ps, s;
    double timeout, um, us, xpos;

    std::string output_file = "rie.txt";

    void starpu(double &p, double &u, const double &mpa) {
        double change, fl, fld, fr, frd, pold, pstart;
        double udiff;
        double tolpre = 1.0e-6;
        int nriter = 20;

        guessp(pstart);
        pold = pstart;
        udiff = ur - ul;

        std::cout << "----------------------------------------\n";
        std::cout << "   Iteration number      Change  \n";
        std::cout << "----------------------------------------\n";

        for (int i = 1; i <= nriter; ++i) {
            prefun(fl, fld, pold, dl, pl, cl);
            prefun(fr, frd, pold, dr, pr, cr);
            p = pold - (fl + fr + udiff) / (fld + frd);
            change = 2.0 * abs((p - pold) / (p + pold));
            std::cout << i << ' ' << change << '\n';
            if (change <= tolpre) {
                goto label;
            }
            if (p < 0.0) {
                p = tolpre;
            }
            pold = p;
        }
        std::cout << "Divergence in Newton-Raphson iteration\n";
    label:
        // Compute velocity in Star Region
        u = 0.5 * (ul + ur + fr - fl);
        
        std::cout << "---------------------------------------\n";
        std::cout << "   Pressure        Velocity\n";
        std::cout << p/mpa << ' ' << u << '\n';
        std::cout << "---------------------------------------\n";
    }

    void guessp(double &pm) {
        double cup, gel, ger, pmax, pmin, ppv, pq;
        double ptl, ptr, qmax, quser, um;
        
        quser = 2.0;
        // Compute guess pressure from PVRS Riemann solver

        cup = 0.25 * (dl + dr) * (cl + cr);
        ppv = 0.5 * (pl + pr) + 0.5 * (ul - ur) * cup;
        ppv = fmax(0.0, ppv);
        pmin = fmin(pl, pr);
        pmax = fmax(pl, pr);
        qmax = pmax / pmin;

        if (qmax <= quser && (pmin <= ppv && ppv <= pmax)) {
            // Select PVRS Riemann solver
            pm = ppv;
        } else {
            if (ppv < pmin) {
                // Select Two-Rarefaction Riemann solver
                pq = pow(pl/pr, g1);
                um = (pq * ul/cl + ur/cr + g4 * (pq - 1.0)) / (pq/cl + 1.0/cr);
                ptl = 1.0 + g7 * (ul - um) / cl;
                ptr = 1.0 + g7 * (um - ur) / cr;
                pm = 0.5 * (pl * pow(ptl, g3) + pr * pow(ptr, g3));
            } else {
                // Select Two-Shock Riemann solver with PVRS as estimate
                gel = sqrt((g5/dl) / (g6 * pl + ppv));
                ger = sqrt((g5/dr) / (g6 * pr + ppv));
                pm = (gel * pl + ger * pr - (ur - ul)) / (gel + ger);
            }
        }
    }

    void prefun(double &f, double &fd, double &p, double &dk, double &pk, double &ck) {
        double ak, bk, prat, qrt;

        if (p <= pk) {
            // Rarefaction wave
            prat = p / pk;
            f = g4 * ck * (pow(prat, g1) - 1.0);
            fd = (1.0/(dk * ck)) * pow(prat, -g2);
        } else {
            // Shock wave
            ak = g5 / dk;
            bk = g6 * pk;
            qrt = sqrt(ak / (bk + p));
            f = (p - pk) * qrt;
            fd = (1.0 - 0.5 * (p - pk) / (bk + p)) * qrt;
        }
    }

    void sample(double &pm, double &um, double &s, double &d, double &u, double &p) {
        double c, cml, cmr, pml, pmr;
        double shl, shr, sl, sr, stl, str;

        if (s <= um) {
            if (pm <= pl) {
                shl = ul - cl;
                if (s <= shl) {
                    d = dl;
                    u = ul;
                    p = pl;
                } else {
                    cml = cl * pow(pm/pl, g1);
                    stl = um - cml;
                    if (s > stl) {
                        d = dl * pow(pm/pl, 1.0/gamma);
                        u = um;
                        p = pm;
                    } else {
                        u = g5 * (cl + g7 * ul + s);
                        c = g5 * (cl + g7 * (ul - s));
                        d = dl * pow(c/cl, g4);
                        p = pl * pow(c/cl, g3);
                    }
                }
            } else {
                // Left shock
                pml = pm / pl;
                sl = ul - cl * sqrt(g2 * pml + g1);
                if (s <= sl) {
                    d = dl;
                    u = ul;
                    p = pl;
                } else {
                    d = dl * (pml + g6) / (pml * g6 + 1.0);
                    u = um;
                    p = pm;
                }
            }
        } else {
            if (pm > pr) {
                // Right shock
                pmr = pm / pr;
                sr = ur + cr * sqrt(g2 * pmr + g1);
                if (s >= sr) {
                    d = dr;
                    u = ur;
                    p = pr;
                } else {
                    d = dr * (pmr + g6) / (pmr * g6 + 1.0);
                    u = um;
                    p = pm;
                }
            } else {
                // Right rarefaction
                shr = ur + cr;
                if (s >= shr) {
                    d = dr;
                    u = ur;
                    p = pr;
                } else {
                    cmr = cr * pow(pm/pr, g1);
                    str = um + cmr;
                    if (s <= str) {
                        d = dr * pow(pm/pr, 1.0/gamma);
                        u = um;
                        p = pm;
                    } else {
                        u = g5 * (-cr + g7 * ur + s);
                        c = g5 * (cr - g7 * (ur - s));
                        d = dr * pow(c/cr, g4);
                        p = pr * pow(c/cr, g3);
                    }
                }
            }
        }
    }

public:
    RiemannSolver(const double &domlen, const double &diaph, const int &cells, const double &gamma, const double &timeout, 
    const double &dl, const double &ul, const double &pl, const double &dr, const double &ur, const double &pr, const double &mpa):
    domlen(domlen), diaph(diaph), cells(cells), gamma(gamma), timeout(timeout),
    dl(dl), ul(ul), pl(pl), dr(dr), ur(ur), pr(pr), mpa(mpa) {
        // Gamma related constants
        g1 = (gamma - 1.0) / (2.0 * gamma);
        g2 = (gamma + 1.0) / (2.0 * gamma);
        g3 = 2.0 * gamma / (gamma - 1.0);
        g4 = 2.0 / (gamma - 1.0);
        g5 = 2.0 / (gamma + 1.0);
        g6 = (gamma - 1.0) / (gamma + 1.0);
        g7 = (gamma - 1.0) / 2.0;
        g8 = gamma - 1.0;

        // Sound speeds
        cl = sqrt(gamma * pl / dl);
        cr = sqrt(gamma * pr / dr);

        // The pressure positivity condition is tested for
        if (g4 * (cl + cr) <= (ur - ul)) {
            // The initial data is such that vaccum is generated
            // Program stopped
            std::cout << "Vacuum is generated by data\n";
            std::cout << "Program stopped\n";
            exit(0);
        }

        // Exact solution for pressure and velocity in star region is found
        starpu(pm, um, mpa);
        dx = domlen / static_cast<double>(cells);

        // Complete solution at time TIMOUT is found
        std::ofstream out;
        out.open(output_file);

        for (int i = 1; i <= cells; ++i) {
            xpos = (static_cast<double>(i) - 0.5) * dx;
            s = (xpos - diaph) / timeout;

            // Solution at point (X, T) = (XPOS - DIAPH, TIMEOUT) is found
            sample(pm, um, s, ds, us, ps);
            // Exact solution profiles are written to rie.txt
            out << xpos << ' ' << ds << ' ' << us << ' ' << ps/mpa << ' ' << ps/ds/g8/mpa << '\n';
        }
        out.close();
    }
};

int main(int argc, char** argv) {
    /*
    if (argc == 1) {
        std::cout << "No input file!\n";
        std::cout << "Example: g++ riemann.cpp input_file\n";
        return 0;
    }
    */
    int cells;
    double domlen, diaph, gamma, timeout;
    double dl, ul, pl, dr, ur, pr, mpa;

    std::string input_file = "exact.ini";
    std::ifstream input;
    input.open(input_file);
    input >> domlen >> diaph >> cells >> gamma >> timeout;
    input >> dl >> ul >> pl >> dr >> ur >> pr >> mpa;
    input.close();

    RiemannSolver rs(domlen, diaph, cells, gamma, timeout, dl, ul, pl, dr, ur, pr, mpa);

    return 0;
} 
