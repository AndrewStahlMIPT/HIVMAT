#ifndef HIVMAT_GAYSS_HPP
#define HIVMAT_GAYSS_HPP

#include <algorithm>

#include <stdexcept>

#include <array>

#include <cstdint>

#include <cmath>

#include <utility>

template < typename RealType, std::uint64_t N >
using Matrix = std::array < std::array < RealType, N > , N > ;

namespace solvers {
    namespace diny {
        template < typename RealType, std::uint32_t N >
        std::pair < Matrix < RealType, N > , std::array < RealType, N >> partialPivot(Matrix < RealType, N > matrix, std::array < RealType, N > b) {
            for (int i = 0; i < N; i++) {
                int pivot_row = i;
                for (int j = i + 1; j < N; j++) {
                    if (std::abs(matrix[j][i]) > std::abs(matrix[pivot_row][i])) {
                        pivot_row = j;
                    }
                }
                if (pivot_row != i) {
                    for (int j = i; j < N; j++) {
                        std::swap(matrix[i][j], matrix[pivot_row][j]);
                    }
                    std::swap(b[i], b[pivot_row]);
                }
                for (int j = i + 1; j < N; j++) {
                    double factor = matrix[j][i] / matrix[i][i];
                    for (int k = i; k < N; k++) {
                        matrix[j][k] -= factor * matrix[i][k];
                    }
                    b[j] -= factor * b[i];
                }
            }
            return std::make_pair(matrix, b);
        }

        template < typename RealType, std::uint32_t N >
        std::array < RealType, N > backSubstitute(const Matrix < RealType, N > & matrix,
                                                  const std::array < RealType, N > & b) {
            std::array < RealType, N > answer;
            for (int i = static_cast < int > (N) - 1; i >= 0; i--) {
                double sum = 0;
                for (std::uint32_t j = i + 1; j < N; j++) {
                    if (j == N) {
                        sum += b[i] * answer[j];
                    } else {
                        sum += matrix[i][j] * answer[j];
                    }
                }
                answer[i] = (b[i] - sum) / matrix[i][i];
            }
            return answer;
        }

    }
    template < typename RealType, std::uint32_t N > [
    [maybe_unused]
    ] std::array < RealType, N > Gauss(const Matrix < RealType, N > & matrix,
                                       const std::array < RealType, N > & b) {
        std::pair < Matrix < RealType, N > , std::array < RealType, N >> temp_pair = diny::partialPivot < RealType, N > (matrix, b);
        return diny::backSubstitute < RealType, N > (temp_pair.first, temp_pair.second);
    }
}
#endif //HIVMAT_GAYSS_HPP