// Utils.hpp
#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <stdexcept>

namespace Utils {
    /**
     * @brief Solves a tridiagonal system of equations using the Thomas Algorithm.
     *
     * @param a Lower diagonal coefficients (size n-1).
     * @param b Main diagonal coefficients (size n).
     * @param c Upper diagonal coefficients (size n-1).
     * @param d Right-hand side vector (size n).
     * @return std::vector<double> Solution vector (size n).
     *
     * @throws std::invalid_argument if vector sizes are incompatible.
     * @throws std::runtime_error if division by zero occurs during the algorithm.
     */
    std::vector<double> thomasAlgorithm(
        const std::vector<double>& a, 
        const std::vector<double>& b, 
        const std::vector<double>& c, 
        const std::vector<double>& d
    );
}

#endif // UTILS_HPP
