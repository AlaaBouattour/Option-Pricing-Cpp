// EuropeanOption.hpp
#ifndef EUROPEANOPTION_HPP
#define EUROPEANOPTION_HPP

#include "Option.hpp"

class EuropeanOption : public Option {
public:
    // Constructors
    EuropeanOption(const std::string& t, const std::string& ex, double mat,
                  double k, const std::string& calcDate, double s,
                  double vol, int m_steps, int n_steps)
        : Option(t, ex, mat, k, calcDate, s, vol, m_steps, n_steps) {}

    // Override price method
    double price(const Market& market) const override;

    // Override clone methods
    EuropeanOption* cloneWithPerturbedSpot(double perturbedSpot) const override;
    EuropeanOption* cloneWithPerturbedVolatility(double perturbedVolatility) const override;
    EuropeanOption* cloneWithPerturbedStrike(double perturbedStrike) const override;
    EuropeanOption* cloneWithPerturbedMaturity(double perturbedMaturity) const override;

    // Override display method
    void display() const override;

    // Greek calculation methods
    double delta(const Market& market, double dS) const;
    double gamma(const Market& market, double dS) const;
    double rho(const Market& market, double dr) const;
    double vega(const Market& market, double dSigma) const;
    double theta(const Market& market, double dT) const;

    // Destructor
    ~EuropeanOption() override {}
};

#endif // EUROPEANOPTION_HPP
