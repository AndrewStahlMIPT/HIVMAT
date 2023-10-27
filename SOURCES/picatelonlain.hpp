#ifndef HIVMAT_PICATELONLAIN_HPP
#define HIVMAT_PICATELONLAIN_HPP

#include <iostream>

#include <iomanip>

#include <cmath>

#include <array>

#include <gtest/gtest.h>

#include <fstream>

template < typename xType, typename yType, unsigned int n >
void data_write([[maybe_unused]] const std::array < xType, n > & xline, [[maybe_unused]]
const std::array < yType, n > & yline,[[maybe_unused]] const std::string&  directory) {

    std::ofstream out;
    out.open(directory + "_Xdata.csv");
    if (out.is_open()) {
        for (auto xelem: xline) {
            out << xelem << std::endl;
        }
    }
    out.close();

    out.open(directory + "_Ydata.csv");
    if (out.is_open()) {
        for (auto yelem: yline) {
            out << yelem << std::endl;

        }
    }
    out.close();

    std::cout << "File has been written" << std::endl;
}


#endif //HIVMAT_PICATELONLAIN_HPP
