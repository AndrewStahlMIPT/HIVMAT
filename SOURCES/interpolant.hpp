#ifndef HIVMAT_INTERPOLANT_HPP
#define HIVMAT_INTERPOLANT_HPP

#include <array>

#include <cmath>

#include <fstream>

#include <iostream>

#include <gtest/gtest.h>

#include<iomanip>

template < typename xType, typename yType, std::size_t N >
class NewtonInterpolator {
private: std::array < xType,
            N > points;
    std::array < yType,
            N > values;

public:
    [[maybe_unused]] NewtonInterpolator(const std::array < xType, N > & points,
                           const std::array < yType, N > & values) noexcept;

    yType interpolate(const xType & x) const noexcept;

};

template < typename xType, typename yType, std::size_t N >
[[maybe_unused]] NewtonInterpolator < xType, yType, N > ::NewtonInterpolator(const std::array < xType, N > & points,
                                                            const std::array < yType, N > & values) noexcept {
    this -> points = points;
    this -> values = values;

    for (std::size_t j = 0; j < N - 1; j++) {
        yType a = this -> values[0];

        for (std::size_t i = 0; i < N - 1 - j; i++) {
            this -> values[i] = (this -> values[i + 1] - this -> values[i]) / (this -> points[i + 1 + j] - this -> points[i]);
        }

        this -> values[N - 1 - j] = a;
    }

}

template < typename xType, typename yType, std::size_t N >
yType NewtonInterpolator < xType, yType, N > ::interpolate(const xType & x) const noexcept {
    yType result = values[0];

    for (std::size_t j = 0; j < N - 1; j++) {
        result *= (x - points[N - 2 - j]);
        result += values[j + 1];
    }

    return result;
}

template < typename xType, unsigned int N >
std::array < xType, N > make_uniform_grid(const xType & a,
                                          const xType & b) {
    xType step = (b - a) / (N - 1);

    std::array < xType, N > grid;

    for (unsigned int i = 0; i < N; i++) {
        grid[i] = a + i * step;
    }

    return grid;
}

template < typename xType, unsigned int N >
std::array < xType, N > make_chebyshev_grid(const xType & a,
                                            const xType & b) {

    std::array < xType, N > grid;

    grid[0] = cos(M_PI / (2 * N));
    xType sin_beta = sin(M_PI / (2 * N));
    xType cos_alpha = cos(M_PI / N);
    xType sin_alpha = sin(M_PI / N);

    for (unsigned int i = 1; i < N; i++) {
        grid[i] = grid[i - 1] * cos_alpha - sin_beta * sin_alpha;
        sin_beta = sin_beta * cos_alpha + grid[i - 1] * sin_alpha;
    }

    for (std::size_t i = 0; i < N; i++) {
        grid[i] = (b + a) / 2 + (b - a) / 2 * grid[i];
    }

    return grid;
}

#endif