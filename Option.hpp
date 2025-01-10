// Option.hpp
#ifndef OPTION_HPP
#define OPTION_HPP

#include <string>

class Market; // Forward declaration

class Option {
public:
    // Members
    std::string type;             // "Call" or "Put"
    std::string exerciseType;     // "européen" or "américain"
    double maturity;              // in years
    double strike;
    std::string calculationDate;
    double spotPrice;
    double volatility;
    int M;                        // Asset price steps
    int N;                        // Time steps

    // Constructor
    Option(const std::string& t, const std::string& ex, double mat,
           double k, const std::string& calcDate, double s,
           double vol, int m_steps, int n_steps)
        : type(t), exerciseType(ex), maturity(mat), strike(k),
          calculationDate(calcDate), spotPrice(s),
          volatility(vol), M(m_steps), N(n_steps) {}

    // Pure virtual methods
    virtual double price(const Market& market) const = 0;

    // Clone methods for Greeks
    virtual Option* cloneWithPerturbedSpot(double perturbedSpot) const = 0;
    virtual Option* cloneWithPerturbedVolatility(double perturbedVolatility) const = 0;
    virtual Option* cloneWithPerturbedStrike(double perturbedStrike) const = 0;
    virtual Option* cloneWithPerturbedMaturity(double perturbedMaturity) const = 0;

    // Pure virtual display method
    virtual void display() const = 0;

    // Virtual destructor
    virtual ~Option() {}
};

#endif // OPTION_HPP
