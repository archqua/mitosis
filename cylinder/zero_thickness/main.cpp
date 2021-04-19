#include "parameters.h"
#include <iostream>
#include <armadillo>

using namespace arma;

vec4(*f)(vec4) = f0;
mat44(*jac)(vec4) = jac0;
void(*initialize)() = initialize0;

mat44 A(vec4 u) {
    return eye(4, 4) / 12 - D * tau / (h*h) - tau / 12 * jac(u); 
}

mat44 B(vec4 u) {
    return -(eye(4, 4) * 5 / 6 + 2 * D * tau / (h*h) - tau * 5 / 6 * jac(u));
}

mat44 C(vec4 u) {
    return eye(4, 4) / 12 - D * tau / (h*h) - tau / 12 * jac(u); 
}

vec4 phi(vec4 u_prev, vec4 u_cur, vec4 u_next) {
    vec4 u_part = u_prev / 12 + u_cur * 5 /6 + u_next / 12; 
    vec4 f_part = tau * (f(u_prev) / 12 + f(u_cur) * 5 / 6 + f(u_next) / 12);
    vec4 jac_part = - tau * (jac(u_prev) * u_prev / 12 + jac(u_cur) * u_cur * 5 / 6 + jac(u_next) * u_next / 12);
    return u_part + f_part + jac_part;
} 

void step(){  
    vec4 new_conc[L];
    mat44 p[L];
    vec4 q[L];

    // A_first * u[0] + B_first * u[1] == phi_first
    mat44 A_first = eye(4,4) * 5 / 6 + 2 * tau / (h*h) * D - tau * jac(conc[0]) * 5 / 6;
    mat44 B_first = eye(4,4) / 6     - 2 * tau / (h*h) * D - tau * jac(conc[1]) / 6;
    vec4 phi_first = conc[0] * 5 / 6 + conc[1] / 6 + tau * (f(conc[0]) * 5 / 6 + f(conc[1]) / 6	- jac(conc[0]) * conc[0] * 5 / 6 - jac(conc[1]) * conc[1] / 6);
    p[0] = - inv(A_first) * B_first;
    q[0] = inv(A_first) * phi_first;

    for (int i = 1; i < L-1; i++) {
        p[i] = - inv(A(conc[i-1]) * p[i-1] - B(conc[i])) * C(conc[i+1]);
        q[i] = inv(A(conc[i-1]) * p[i-1] - B(conc[i])) * (phi(conc[i-1], conc[i], conc[i+1]) - A(conc[i-1]) * q[i-1]);
    }

    // A_last * u[L-2] + B_last * u[L-1] == phi_last
    mat44 A_last = eye(4, 4) / 6     - 2 * tau * D / (h*h) - tau * jac(conc[L-2]) / 6;
    mat44 B_last = eye(4, 4) * 5 / 6 + 2 * tau * D / (h*h) - tau * jac(conc[L-1]) * 5 / 6;
    vec4 phi_last = conc[L-2] / 6 + conc[L-1] * 5 / 6 + tau * (f(conc[L-2]) / 6 + f(conc[L-1]) * 5 / 6 - jac(conc[L-2]) * conc[L-2] / 6 - jac(conc[L-1]) * conc[L-1] * 5 / 6);
    new_conc[L-1] = inv(A_last * p[L-2] + B_last) * (phi_last - A_last * q[L-2]);
    for (int i = L - 2; i > -1; i--){
        new_conc[i] = p[i] * new_conc[i+1] + q[i];
    }

    std::copy(new_conc, new_conc+L, conc);
}

int main(int argc, char *argv[]) {
    
    bool use_test = false;
    char n_test = 0;
    int c;
    while ((c = getopt(argc, argv, "t:")) != -1) {
        switch (c) {
            case 't':
                use_test = true;
                n_test = atoi(optarg);
                break;
            case '?':
                break;
            default:
                break;
        }
    }
    if (use_test) {
        f = f_test[n_test];
        jac = jac_test[n_test]; 
        initialize = initialize_test[n_test];
    }

    initialize();
    step();    
    
    for (auto &v : conc) {
        std::cout << v << std::endl;
    }

    return 0;
}
