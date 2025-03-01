// AmericanOption.hpp
#ifndef AMERICANOPTION_HPP
#define AMERICANOPTION_HPP

#include "Option.hpp"

class AmericanOption : public Option {
public:
    // Constructors
    AmericanOption(const std::string& t, const std::string& ex, double mat,
                  double k, const std::string& calcDate, double s,
                  double vol, int m_steps, int n_steps)
        : Option(t, ex, mat, k, calcDate, s, vol, m_steps, n_steps) {}

    // Override price method
    double price(const Market& market) const override;

    // Override clone methods
    AmericanOption* cloneWithPerturbedSpot(double perturbedSpot) const override;
    AmericanOption* cloneWithPerturbedVolatility(double perturbedVolatility) const override;
    AmericanOption* cloneWithPerturbedStrike(double perturbedStrike) const override;
    AmericanOption* cloneWithPerturbedMaturity(double perturbedMaturity) const override;

    // Override display method
    void display() const override;

    // Greek calculation methods
    double delta(const Market& market, double dS) const;
    double gamma(const Market& market, double dS) const;
    double rho(const Market& market, double dr) const;
    double vega(const Market& market, double dSigma) const;
    double theta(const Market& market, double dT) const;

    // Destructor
    ~AmericanOption() override {}
};

#endif // AMERICANOPTION_HPP
