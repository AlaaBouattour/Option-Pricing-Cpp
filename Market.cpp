// Market.cpp
#include "Market.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

void Market::loadRates(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open rates file.");
    }

    std::string line;
    // Skip header
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double time, rate;
        std::string temp;
        std::getline(ss, temp, ','); // Time
        time = std::stod(temp);
        std::getline(ss, temp, ','); // Rate
        rate = std::stod(temp);
        rates[time] = rate;
    }

    file.close();
}

double Market::getRate(double t) const {
    if (rates.empty()) {
        throw std::runtime_error("No rates loaded.");
    }

    // If exact time exists
    auto it = rates.find(t);
    if (it != rates.end()) {
        return it->second;
    }

    // Otherwise, perform linear interpolation
    auto upper = rates.upper_bound(t);
    if (upper == rates.begin()) {
        return upper->second;
    }
    if (upper == rates.end()) {
        return std::prev(upper)->second;
    }

    auto lower = std::prev(upper);
    double t1 = lower->first;
    double r1 = lower->second;
    double t2 = upper->first;
    double r2 = upper->second;

    // Linear interpolation
    double interpolatedRate = r1 + (r2 - r1) * (t - t1) / (t2 - t1);
    return interpolatedRate;
}

void Market::setRate(double t, double newRate) {
    rates[t] = newRate;
}

void Market::displayRates() const {
    for (const auto& [time, rate] : rates) {
        std::cout << "Time: " << time << " years, Rate: " << rate << "\n";
    }
}
