#ifndef COMPUTEYASSPARAMETERS_H
#define COMPUTEYASSPARAMETERS_H

#include <cstdlib>
#include <cmath>
#include <vector>



inline double dpow(double a, size_t i_expon)
{
    double c_2px = a;
    double result = 1;

    while (i_expon) {
      if (i_expon & (size_t) 1) {
        result *= c_2px;
      }
      c_2px *= c_2px;
      i_expon >>= 1;
    }
    return result;
}


inline size_t computeRho(double mutationRate, size_t weight, double alpha) {
    auto k = weight;
    double p = 1. - mutationRate;

    double p_k = dpow(p, k);
    double qp_k = (1 - p) * p_k;
    size_t x = 0;
    double Sum_0_xk1 = 0.00;
    double Sum = 0;

    auto last_k_prob = std::vector<double>(k+1);

    while (Sum < (1 - alpha)) {
        if (x < k) {
            last_k_prob[x % (k + 1)] = 0.00;
        } else if (x == k) {
            last_k_prob[x % (k + 1)] = p_k;
        } else {
            Sum_0_xk1 += last_k_prob[x % (k + 1)];
            last_k_prob[x % (k + 1)] = qp_k * (1 - Sum_0_xk1);
        }

        Sum += last_k_prob[x % (k + 1)];
        x++;
    }
    //FREE(last_k_prob, (k + 1) * sizeof(double));

    return x;
}



inline std::vector<double> randomwalk_probability_of_pos3(double pI, size_t L) {
    size_t i, j;
    auto P = L;
    double a = pI * 0.50;
    double b = 1 - pI;
    auto u = std::vector<double>(2*L+1);
    auto f = std::vector<double>(2*L+1);
    auto t = std::vector<double>(2*L+1);
    std::vector<double> s;

    /* (1) tables inits */

    for (i = 0; i < 2 * L + 1; i++)
        u[i] = 0;
    for (i = 0; i < 2 * L + 1; i++)
        f[i] = 0;

    f[0] = 1.00;
    u[0] = a;
    u[1] = b;
    u[2] = a;

    /* (2) qpow */
    while (P > 0) {
        if (P & 1) {
            /* f = f*u */
            for (i = 0; i < 2 * L + 1; i++) {
                t[i] = 0;
                for (j = 0; j <= i; j++)
                    t[i] += u[i - j] * f[j];
            }
            s = t;
            t = f;
            f = s;
        }
        /* u = u * u; */
        for (i = 0; i < 2 * L + 1; i++) {
            t[i] = 0;
            for (j = 0; j <= i; j++)
                t[i] += u[i - j] * u[j];
        }
        s = t;
        t = u;
        u = s;
        P >>= 1;
    }

    return f;

}



inline size_t computeDelta(double indelRate, size_t rho, double alpha) {
    auto pI = indelRate;
    auto L = rho;

    auto RDW_Bound = randomwalk_probability_of_pos3(pI, L);
    double Sum = RDW_Bound[L];
    size_t bound = 0;
    while (Sum < (1. - alpha) && bound < L) {
        bound++;
        Sum += RDW_Bound[L - bound] + RDW_Bound[L + bound];
    }

    return bound;
}

#endif // COMPUTEYASSPARAMETERS_H
