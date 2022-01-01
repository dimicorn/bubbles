#include <string>
#ifndef CLASS_H
#define CLASS_H

class Bubble {
    private:
        double gamma_, k_rho_, n_int_;
        double eta_ = (2 + n_int_) / (5 - k_rho_);

    public:
        explicit Bubble(double gamma, double k_rho, double n_int);
        double CurveValue(double lambda_c);
        double LambdaApprox();
        double GradVel1();
        double GradVel2();
        double Delta(double lambda_c, double lambda_c_approx);
        double R_sw(double vel, double numb, double k);
        int Norm1();
        int Norm2();
        int Values();
        int Q_p();
        std::string Output();
};
#endif
