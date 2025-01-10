// Market.hpp
#ifndef MARKET_HPP
#define MARKET_HPP

#include <map>
#include <string>

class Market {
private:
    std::map<double, double> rates; // Mapping from time to rate

public:
    // Load rates from a CSV file
    void loadRates(const std::string& filename);

    // Get interpolated rate at time t
    double getRate(double t) const;

    // Set rate at a specific time (for Rho calculations)
    void setRate(double t, double newRate);

    // Optional: Display rates
    void displayRates() const;
};

#endif // MARKET_HPP
