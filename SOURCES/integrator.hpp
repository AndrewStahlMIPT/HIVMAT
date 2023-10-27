#ifndef HIVMAT_INTEGRATOR_HPP
#define HIVMAT_INTEGRATOR_HPP

#include <iostream>

#include <random>

#include <cmath>

#include <array>

#include <type_traits>

#include <vector>

const int n = 3;

template < typename RealType >
RealType Function(RealType x) {
    return std::sin(x);
}

std::array < std::array < double, n > , n > get_coefficients(int _pl, std::array < double, n > _xi) {

    std::array < std::array < double, 2 > , n > coefficients {};
    for (int i = 0; i < n; i++) {
        if (i == _pl) {
            coefficients[i][0] = std::numeric_limits < double > ::infinity();
            coefficients[i][1] = std::numeric_limits < double > ::infinity();
        } else {
            coefficients[i][0] = 1 / (_xi[_pl] - _xi[i]);
            coefficients[i][1] = -_xi[i] / (_xi[_pl] - _xi[i]);
        }
    }

    std::array < std::array < double, n > , n > filtered_coefficients {};
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (coefficients[i][0] != std::numeric_limits < double > ::infinity()) {
            filtered_coefficients[j][0] = coefficients[i][1];
            filtered_coefficients[j][1] = coefficients[i][0];
            j += 1;
        }
    }
    return filtered_coefficients;
}

std::array < std::array < double, n > , n > get_polynomial_l(const std::array < double, n > _xi) {
    std::array < std::array < double, n > , n > pli {};
    for (int pl = 0; pl < n; pl++) {
        std::array < std::array < double, n > , n > coefficients = get_coefficients(pl, _xi);
        for (int i = 1; i < n - 1; i++) {
            if (i == 1) {
                pli[pl][0] = coefficients[i - 1][0] * coefficients[i][0];
                pli[pl][1] = coefficients[i - 1][1] * coefficients[i][0] + coefficients[i][1] * coefficients[i - 1][0];
                pli[pl][2] = coefficients[i - 1][1] * coefficients[i][1];
            } else {
                std::array < double ,n> clone_pli{};
                for (int val = 0; val < n; val++) {
                    clone_pli[val] = pli[pl][val];
                }
                std::array < double,n > zeros_pli{};
                for (int j = 0; j < n - 1; j++) {
                    double product_1 = clone_pli[j] * coefficients[i][0];
                    double product_2 = clone_pli[j] * coefficients[i][1];
                    zeros_pli[j] += product_1;
                    zeros_pli[j + 1] += product_2;
                }
                for (int val = 0; val < n; val++) {
                    pli[pl][val] = zeros_pli[val];
                }
            }
        }
    }
    return pli;
}

std::array < double, n > getEasyWeights(double minus_one, double plus_one,
                                        const std::array < double, n > uzelki) {
    std::array < double, n > easyWeights {};
    for (int i = 0; i < uzelki.size(); i++) {
        for (int j = 0; j < uzelki.size(); j++) {
            easyWeights[i] += (std::pow(plus_one, j + 1) - std::pow(minus_one, j + 1)) / (j + 1) * get_polynomial_l(uzelki)[i][j];
        }
    }
    return easyWeights;
}

template < typename RealType >
RealType integratorGauss(const RealType a,
                         const RealType b,
                         const std::array < RealType, n > uzelki, RealType h) {

    const std::array < RealType, n > w = getEasyWeights(-1, 1, uzelki);
    RealType I = 0, dI = 0;
    RealType a_i, b_i, arg1, arg2, xi_;

    for (auto i = 0; i < (b-a)/h; ++i) {
        a_i = a + i * h;
        b_i = a + (i + 1) * h;
        arg1 = (b_i - a_i) / 2;
        arg2 = (a_i + b_i) / 2;
        for (auto j = 0; j < n; ++j) {
            xi_ = arg2 + arg1 * uzelki[j];
            dI += w[j] * Function(xi_);
        }
        dI *= arg1;
        I += dI;
        dI = 0;
    }
    return I;
}

#endif //HIVMAT_INTEGRATOR_HPP