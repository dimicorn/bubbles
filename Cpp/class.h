#ifndef CLASS_H
#define CLASS_H

class Bubble {
    private:
        float gamma_, k_rho_, n_int_;
        float eta_ = (2 + n_int_) / (5 - k_rho_);

    public:
        Bubble(float Gamma, float k_rho, float n_int);
        float curve_value(float lambda_c);
        float lambda_approx();
        float dv1();
        float dv2();
        float delta(float lambda_c, float lambda_c_approx);
        float r_sw(float vel, float numb, float k);
        int norm1();
        int norm2();
        int values();
};
#endif
