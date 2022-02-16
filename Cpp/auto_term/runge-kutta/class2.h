#ifndef CLASS2_H
#define CLASS2_H
class Bubble {
    private:
        double gamma__, k_rho__, n_int__;
        double eta__ = (2 + n_int__) / (5 - k_rho__);
        int i_, j_, k_;

    public:
        explicit Bubble(double gamma, double k_rho, double n_int, int i, int j, int l);
        double CurveValue(double lambda_c);
        double LambdaApprox();
};
#endif
