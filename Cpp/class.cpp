#include <cmath>
#include "class.h"
#include "constants.h"

Bubble::Bubble(float Gamma, float k_rho, float n_int):
    gamma_(Gamma), k_rho_(k_rho), n_int_(n_int) {
    // Solving system of ODE's
    // Checking for intersection
    // Other magic happening 
    };

// Value of the curve at lambda_c
float Bubble::curve_value(float lambda_c) {
    return (gamma_ + 1) / 2 * lambda_c;
}

// Approximation using eqn B8a
float Bubble::lambda_approx() {
    float t = gamma_ * gamma_ * gamma_ + 12 * gamma_ * gamma_ + 8 * gamma_ + 1 - 0.5 * (gamma_ + 1) * (3 * gamma_ + 1) * k_rho_ - (gamma_ + 1) * (4 * gamma_ + 1) / eta_;
    float u = 3 * gamma_ * gamma_ * gamma_ + 12 * gamma_ * gamma_ + 7 * gamma_ + 1 - 0.5 * (gamma_ + 1) * (3 * gamma_ + 1) * k_rho_ - (gamma_ + 1) * (4 * gamma_ + 1) / eta_;
    return t / u; 
}

// Gradient of velocity, r = R_s
float Bubble::dv1() {
    float temp = (-(7 * gamma_ + 3) + (gamma_ + 1) * k_rho_ + 3 * (gamma_ + 1) / eta_) / (gamma_ + 1);
    return temp;
}

// Gradient of velocity, r = R_s
float Bubble::dv2() {
    float temp = (-2 * (gamma_ + 1) + k_rho_ + 2 / eta_) / (2 * gamma_ / (gamma_ + 1));
    return temp;
}

// Difference between lambda from the ODE and the approximation using eqn B8a
float Bubble::delta(float lambda_c, float lambda_c_approx) {
    float temp = abs(lambda_c - lambda_c_approx) / lambda_c * 100;
    return temp;
}

// Calculating R_sw
float Bubble::r_sw(float vel, float numb, float k) {
    float f_rho = pow((4 * gamma_ / (gamma_ + 1) * (gamma_ + 1)), 1 / (gamma_ - 1));
    float time = 3600 * 24 * 365 * numb;
    float temp = (gamma_ - 1) / (gamma_ + 1) * vel * f_rho * time * (3 * gamma_ / (3 * (gamma_ - 1) * eta_ + n_int_)) / (k_rho_ * k_rho_ * k_rho_);
    return temp / au;
}

int Bubble::values() {
    return 0;
}

int Bubble::norm1() {
    return 0;
}

int Bubble::norm2() {
    return 0;
}

int Bubble::q_p() {
     return 0;
}
