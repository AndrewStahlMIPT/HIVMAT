#include <gtest/gtest.h>

#include "../SOURCES/num_dif.hpp"

#include "../SOURCES/picatelonlain.hpp"

#include <iostream>

#include <iomanip>

#include <array>

#include <cmath>

TEST(test_simp_1_numdif, N2) {

    const unsigned int N = 2;
    std::array < double, N > points = {
            -1,
            1
    };
    DerivativeCoef < double, N > test = calcDerivativeCoef < double, N > (points);
    std::cout << test.centralCoef << '\n';
    for (unsigned int i = 0; i < N; i++) {
        std::cout << test.otherCoefs[i] << '\n';
    }
}

TEST(test_simp_2_numdif, N2) {

    const unsigned int N = 2;
    std::array < double, N > points = {
            1,
            2
    };
    DerivativeCoef < double, N > test = calcDerivativeCoef < double, N > (points);
    std::cout << test.centralCoef << '\n';
    for (unsigned int i = 0; i < N; i++) {
        std::cout << test.otherCoefs[i] << '\n';
    }
}

TEST(test_numdif_1, N3) {
    const unsigned int x = 1;
    const double f_x = std::exp(x);
    const unsigned int N = 3;
    unsigned int count = 0;
    double error;
    double h = 1e-15;

    std::array < double, 16 > H {};
    std::array < double, 16 > Error {};
    std::array < double, N > points = {};
    std::array < double, N > func {};

    while (h <= 10) {
        for (unsigned int i = 0; i < N; i++) {
            points[i] = h * (i + 1);
            func[i] = std::exp(x + (i + 1) * h);
        }

        DerivativeCoef < double, N > test = calcDerivativeCoef < double, N > (points);
        double proizv = f_x * test.centralCoef;

        for (unsigned int i = 0; i < N; i++) {
            proizv = test.otherCoefs[i] * func[i] + proizv;
        }
        error = std::abs(f_x - proizv);
        Error[count] = error;
        H[count] = h;
        h *= 10;
        count++;

    }
    data_write < double, double, 16 > (H, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/NUM_DIF_N3");
}

TEST(test_numdif_2, N4) {
    const unsigned int x = 1;
    const double f_x = std::exp(x);
    const unsigned int N = 4;
    unsigned int count = 0;
    double error;
    double h = 1e-15;

    std::array < double, 16 > H {};
    std::array < double, 16 > Error {};
    std::array < double, N > points = {};
    std::array < double, N > func {};

    while (h <= 10) {
        for (unsigned int i = 0; i < N; i++) {
            points[i] = h * (i + 1);
            func[i] = std::exp(x + (i + 1) * h);
        }

        DerivativeCoef < double, N > test = calcDerivativeCoef < double, N > (points);
        double proizv = f_x * test.centralCoef;

        for (unsigned int i = 0; i < N; i++) {
            proizv = test.otherCoefs[i] * func[i] + proizv;
        }
        error = std::abs(f_x - proizv);
        Error[count] = error;
        H[count] = h;
        h *= 10;
        count++;

    }
    data_write < double, double, 16 > (H, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/NUM_DIF_N4");
}

TEST(test_numdif_3, N5) {
    const unsigned int x = 1;
    const double f_x = std::exp(x);
    const unsigned int N = 5;
    unsigned int count = 0;
    double error;
    double h = 1e-15;

    std::array < double, 16 > H {};
    std::array < double, 16 > Error {};
    std::array < double, N > points = {};
    std::array < double, N > func {};

    while (h <= 10) {
        for (unsigned int i = 0; i < N; i++) {
            points[i] = h * (i + 1);
            func[i] = std::exp(x + (i + 1) * h);
        }

        DerivativeCoef < double, N > test = calcDerivativeCoef < double, N > (points);
        double proizv = f_x * test.centralCoef;

        for (unsigned int i = 0; i < N; i++) {
            proizv = test.otherCoefs[i] * func[i] + proizv;
        }
        error = std::abs(f_x - proizv);
        Error[count] = error;
        H[count] = h;
        h *= 10;
        count++;

    }
    data_write < double, double, 16 > (H, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/NUM_DIF_N5");
}