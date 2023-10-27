#include <gtest/gtest.h>

#include "../SOURCES/interpolant.hpp"

#include "../SOURCES/picatelonlain.hpp"

#include <cmath>

#include <array>

TEST(test_1, N3) {

    double a = 0, b = 1.0 / 16;
    const unsigned int samp_freq = 1001;
    const unsigned int N = 3;
    unsigned int count = 0;
    std::array < double, 6 > B{};
    std::array < double, 6 > Error{};
    while (b <= 2) {
        std::array < double, N > x = make_uniform_grid < double, N > (a, b);
        std::array < double, N > y = x;
        std::array < double, samp_freq > yfunc {};
        std::array < double, samp_freq > xfunc = make_uniform_grid < double, samp_freq > (a, b);
        for (unsigned long i = 0; i < yfunc.size(); i++) {
            yfunc[i] = exp(xfunc[i]);
        }
        for (double & i : y) {
            i = exp(i);
        }
        NewtonInterpolator < double, double, N > interpolator(x, y);
        double error = -1;
        for (unsigned long i = 0; i < xfunc.size(); i++) {
            double result = interpolator.interpolate(xfunc[i]);
            error = std::max(abs(result - yfunc[i]), error);
            Error[count] = error;
        }
        B[count] = b;
        b *= 2;
        count++;
    }
    data_write<double, double, 6>(B, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/N3");
}

TEST(test_2, N4) {

    double a = 0, b = 1.0 / 16;
    const unsigned int samp_freq = 1001;
    const unsigned int N = 4;
    unsigned int count = 0;
    std::array < double, 6 > B{};
    std::array < double, 6 > Error{};
    while (b <= 2) {
        std::array < double, N > x = make_uniform_grid < double, N > (a, b);
        std::array < double, N > y = x;
        std::array < double, samp_freq > yfunc {};
        std::array < double, samp_freq > xfunc = make_uniform_grid < double, samp_freq > (a, b);
        for (unsigned long i = 0; i < yfunc.size(); i++) {
            yfunc[i] = exp(xfunc[i]);
        }
        for (double & i : y) {
            i = exp(i);
        }
        NewtonInterpolator < double, double, N > interpolator(x, y);
        double error = -1;
        for (unsigned long i = 0; i < xfunc.size(); i++) {
            double result = interpolator.interpolate(xfunc[i]);
            error= std::max(abs(result - yfunc[i]), error);
            Error[count] = error;
        }
        B[count] = b;
        b *= 2;
        count++;
    }
    data_write<double, double, 6>(B, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/N4");
}

TEST(test_3, N5) {

    double a = 0, b = 1.0 / 16;
    const unsigned int samp_freq = 1001;
    const unsigned int N = 5;
    unsigned int count = 0;
    std::array < double, 6 > B{};
    std::array < double, 6 > Error{};
    while (b <= 2) {
        std::array < double, N > x = make_uniform_grid < double, N > (a, b);
        std::array < double, N > y = x;
        std::array < double, samp_freq > yfunc {};
        std::array < double, samp_freq > xfunc = make_uniform_grid < double, samp_freq > (a, b);
        for (unsigned long i = 0; i < yfunc.size(); i++) {
            yfunc[i] = exp(xfunc[i]);
        }
        for (double & i : y) {
            i = exp(i);
        }
        NewtonInterpolator < double, double, N > interpolator(x, y);
        double error = -1;
        for (unsigned long i = 0; i < xfunc.size(); i++) {
            double result = interpolator.interpolate(xfunc[i]);
            error = std::max(abs(result - yfunc[i]), error);
            Error[count] = error;
        }
        B[count] = b;
        b *= 2;
        count++;
    }
    data_write<double, double, 6>(B, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/N5");
}

TEST(test_4, N3) {

    double a = 0, b = 1.0 / 16;
    const unsigned int samp_freq = 1001;
    const unsigned int N = 3;
    unsigned int count = 0;
    std::array < double, 6 > B{};
    std::array < double, 6 > Error{};
    while (b <= 2) {
        std::array < double, N > x = make_chebyshev_grid < double, N > (a, b);
        std::array < double, N > y = x;
        std::array < double, samp_freq > yfunc {};
        std::array < double, samp_freq > xfunc = make_chebyshev_grid < double, samp_freq > (a, b);
        for (unsigned long i = 0; i < yfunc.size(); i++) {
            yfunc[i] = exp(xfunc[i]);
        }
        for (double & i : y) {
            i = exp(i);
        }
        NewtonInterpolator < double, double, N > interpolator(x, y);
        double error = -1;
        for (unsigned long i = 0; i < xfunc.size(); i++) {
            double result = interpolator.interpolate(xfunc[i]);
            error = std::max(abs(result - yfunc[i]), error);
            Error[count] = error;
        }
        B[count] = b;
        b *= 2;
        count++;
    }
    data_write<double, double, 6>(B, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/N3CH");
}

TEST(test_5, N4) {

    double a = 0, b = 1.0 / 16;
    const unsigned int samp_freq = 1001;
    const unsigned int N = 4;
    unsigned int count = 0;
    std::array < double, 6 > B{};
    std::array < double, 6 > Error{};
    while (b <= 2) {
        std::array < double, N > x = make_chebyshev_grid < double, N > (a, b);
        std::array < double, N > y = x;
        std::array < double, samp_freq > yfunc {};
        std::array < double, samp_freq > xfunc = make_chebyshev_grid < double, samp_freq > (a, b);
        for (unsigned long i = 0; i < yfunc.size(); i++) {
            yfunc[i] = exp(xfunc[i]);
        }
        for (double & i : y) {
            i = exp(i);
        }
        NewtonInterpolator < double, double, N > interpolator(x, y);
        double error = -1;
        for (unsigned long i = 0; i < xfunc.size(); i++) {
            double result = interpolator.interpolate(xfunc[i]);
            error = std::max(abs(result - yfunc[i]), error);
            Error[count] = error;
        }
        B[count] = b;
        b *= 2;
        count++;
    }
    data_write<double, double, 6>(B, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/N4CH");
}

TEST(test_6, N5) {

    double a = 0, b = 1.0 / 16;
    const unsigned int samp_freq = 1001;
    const unsigned int N = 5;
    unsigned int count = 0;
    std::array < double, 6 > B{};
    std::array < double, 6 > Error{};
    while (b <= 2) {
        std::array < double, N > x = make_chebyshev_grid < double, N > (a, b);
        std::array < double, N > y = x;
        std::array < double, samp_freq > yfunc {};
        std::array < double, samp_freq > xfunc = make_chebyshev_grid < double, samp_freq > (a, b);
        for (unsigned long i = 0; i < yfunc.size(); i++) {
            yfunc[i] = exp(xfunc[i]);
        }
        for (double & i : y) {
            i = exp(i);
        }
        NewtonInterpolator < double, double, N > interpolator(x, y);
        double error = -1;
        for (unsigned long i = 0; i < xfunc.size(); i++) {
            double result = interpolator.interpolate(xfunc[i]);
            error = std::max(abs(result - yfunc[i]), error);
            Error[count] = error;
        }
        B[count] = b;
        b *= 2;
        count++;
    }
    data_write<double, double, 6>(B, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/N5CH");
}