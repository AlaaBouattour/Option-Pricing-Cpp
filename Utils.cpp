// Utils.cpp
#include "Utils.hpp"
#include <vector>
#include <stdexcept>

namespace Utils {
    std::vector<double> thomasAlgorithm(
        const std::vector<double>& a, 
        const std::vector<double>& b, 
        const std::vector<double>& c, 
        const std::vector<double>& d
    ) {
        int n = b.size();
        if (a.size() != n - 1 || c.size() != n - 1 || d.size() != n) {
            throw std::invalid_argument("Invalid vector sizes for Thomas Algorithm.");
        }

        std::vector<double> c_prime(n - 1, 0.0);
        std::vector<double> d_prime(n, 0.0);
        std::vector<double> x(n, 0.0);

        // Forward sweep
        c_prime[0] = c[0] / b[0];
        d_prime[0] = d[0] / b[0];

        for(int i = 1; i < n -1; ++i){
            double m = b[i] - a[i-1] * c_prime[i-1];
            if(m == 0.0){
                throw std::runtime_error("Division by zero encountered in Thomas Algorithm.");
            }
            c_prime[i] = c[i] / m;
            d_prime[i] = (d[i] - a[i-1] * d_prime[i-1]) / m;
        }

        // Last element of d_prime
        double m_last = b[n-1] - a[n-2] * c_prime[n-2];
        if(m_last == 0.0){
            throw std::runtime_error("Division by zero encountered in Thomas Algorithm at last element.");
        }
        d_prime[n-1] = (d[n-1] - a[n-2] * d_prime[n-2]) / m_last;

        // Back substitution
        x[n-1] = d_prime[n-1];
        for(int i = n -2; i >=0; --i){
            x[i] = d_prime[i] - c_prime[i] * x[i+1];
        }

        return x;
    }
}
