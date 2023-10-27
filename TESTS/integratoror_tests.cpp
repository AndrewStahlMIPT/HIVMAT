
#include "gtest/gtest.h"

#include "../SOURCES/integrator.hpp"

#include "../SOURCES/picatelonlain.hpp"

TEST(test_integrator_1, N3) {

    const double N = 10.;
    const double a = 0.,
                 b = 10.;
    double h = (b-a)/N;
    double error;
    unsigned int count = 0;
    double  Int_ist = 1.839071529076452;
    std::array < double, 6 > H {};
    std::array < double, 6 > Error {};
    std::array < double, n > uzelki = {
            -std::sqrt(3.0 / 5.0),
            0,
            std::sqrt(3.0 / 5.0)
    };
    while (h >= 0.00001) {
        error = std::abs(Int_ist -integratorGauss(a, b, uzelki, h));
        Error[count] = error;
        H[count] = h;
        h /= 10;
        count++;
    }
    data_write < double, double, 6 > (H, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/NUM_INTEGRATE_N3_1");
}

TEST(test_integrator_2, N3) {

    const double N = 10.;
    const double a = 0.,
            b = 10.;
    double h = (b-a)/N;
    double error;
    unsigned int count = 0;
    double  Int_ist = 1.839071529076452;
    std::array < double, 6 > H {};
    std::array < double, 6 > Error {};
    std::array < double, n > uzelki = {
            -std::sqrt(4.0 / 5.0),
            0,
            std::sqrt(4.0 / 5.0)
    };
    while (h >= 0.00001) {
        error = std::abs(Int_ist -integratorGauss(a, b, uzelki, h));
        Error[count] = error;
        H[count] = h;
        h /= 10;
        count++;
    }
    data_write < double, double, 6 > (H, Error, "/Users/andrewstahl/CLionProjects/HIVMAT/GRAPHS/NUM_INTEGRATE_N3_2");
}