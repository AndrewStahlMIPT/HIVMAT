#ifndef HIVMAT_NUM_DIF_HPP
#define HIVMAT_NUM_DIF_HPP

#include <array>

#include <cmath>

#include <fstream>

#include <iostream>

#include <gtest/gtest.h>

#include <iomanip>

#include "../SOURCES/GAYSS.hpp"

template < typename RealType, unsigned int N >
struct DerivativeCoef {
    RealType centralCoef;
    std::array < RealType, N > otherCoefs;
};

template < typename RealType, unsigned int N >
DerivativeCoef < RealType, N > calcDerivativeCoef(const std::array < RealType, N > & points) noexcept {
    std::array < std::array < RealType, N + 1 > , N + 1 > matrix;
    std::array < RealType, N + 1 > b {};
    b[1] = 1;
    for (unsigned int i = 0; i < N + 1; i++) {
        matrix[0][i] = 1;
    }
    for (unsigned int i = 0; i < N; i++) {
        matrix[i + 1][0] = 0;
        for (unsigned int k = 0; k < N; k++) {
            matrix[i + 1][k + 1] = matrix[i][k + 1] * points[k];
            matrix[i + 1][k + 1] /= (i + 1);
        }
    }
    std::array < RealType, N + 1 > b_final = solvers::Gauss < RealType, N + 1 > (matrix, b);
    DerivativeCoef < RealType, N > iskomoe;
    iskomoe.centralCoef = b_final[0];
    std::array < RealType, N > haha;
    for (unsigned int i = 0; i < N; i++) {
        haha[i] = b_final[i + 1];
    }
    iskomoe.otherCoefs = haha;
    return iskomoe;
}

#endif //HIVMAT_NUM_DIF_HPP