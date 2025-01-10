// EuropeanOption.cpp
#include "EuropeanOption.hpp"
#include "Market.hpp"   
#include "Utils.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>

// On utilise les clones pour plus d'efficacité sur la mémoire
EuropeanOption* EuropeanOption::cloneWithPerturbedSpot(double perturbedSpot) const {
    return new EuropeanOption(type, exerciseType, maturity, strike,
                              calculationDate, perturbedSpot, volatility, M, N);
}

EuropeanOption* EuropeanOption::cloneWithPerturbedVolatility(double perturbedVolatility) const {
    return new EuropeanOption(type, exerciseType, maturity, strike,
                              calculationDate, spotPrice, perturbedVolatility, M, N);
}

EuropeanOption* EuropeanOption::cloneWithPerturbedStrike(double perturbedStrike) const {
    return new EuropeanOption(type, exerciseType, maturity, perturbedStrike,
                              calculationDate, spotPrice, volatility, M, N);
}

EuropeanOption* EuropeanOption::cloneWithPerturbedMaturity(double perturbedMaturity) const {
    return new EuropeanOption(type, exerciseType, perturbedMaturity, strike,
                              calculationDate, spotPrice, volatility, M, N);
}

// Affichage
void EuropeanOption::display() const {
    std::cout << "European Option Details:\n";
    std::cout << "Type: " << type << "\n";
    std::cout << "Exercise Type: " << exerciseType << "\n";
    std::cout << "Maturity: " << maturity << " years\n";
    std::cout << "Strike Price: " << strike << "\n";
    std::cout << "Calculation Date: " << calculationDate << "\n";
    std::cout << "Spot Price: " << spotPrice << "\n";
    std::cout << "Volatility: " << volatility << "\n";
    std::cout << "Asset Price Steps (M): " << M << "\n";
    std::cout << "Time Steps (N): " << N << "\n";
}

// Finite Difference Method
double EuropeanOption::price(const Market& market) const {
    // Retrieve the risk-free rate at current time
    double r = market.getRate(maturity);

    // Asset price range
    double S_max = spotPrice * 3.0; 
    double dS = S_max / M;

    // Time step
    double dt = maturity / N;

    // Initialize asset price grid
    std::vector<double> S(M + 1);
    for(int i = 0; i <= M; ++i){
        S[i] = i * dS;
    }

    // Initialize option values at maturity
    std::vector<double> V_prev(M + 1, 0.0);
    for(int i = 0; i <= M; ++i){
        if(type == "Call"){
            V_prev[i] = std::max(S[i] - strike, 0.0);
        }
        else if(type == "Put"){
            V_prev[i] = std::max(strike - S[i], 0.0);
        }
        else{
            throw std::runtime_error("Invalid option type.");
        }
    }

    // Number of interior nodes
    int n = M - 1;

    // Initialize coefficient vectors
    std::vector<double> a_coeff(n - 1, 0.0); // Lower diagonal (a1 to a_{n-1})
    std::vector<double> b_coeff(n, 0.0);     // Main diagonal (b1 to bn)
    std::vector<double> c_coeff(n - 1, 0.0); // Upper diagonal (c1 to c_{n-1})
    std::vector<double> d_coeff(n, 0.0);     // Right-hand side (d1 to dn)

    // Initialize V_new
    std::vector<double> V_new(M + 1, 0.0);

    // Variables to store alpha for first and gamma for last node
    double alpha_first = 0.0;
    double gamma_last = 0.0;

    // Time-stepping loop
    for(int step = N -1; step >=0; --step){
        // Populate the tridiagonal matrix coefficients
        for(int j = 1; j <=n; ++j){ // j =1 to n
            double Si = S[j];
            double alpha = 0.25 * dt * (volatility * volatility * Si * Si / (dS * dS) - r * Si / dS);
            double beta  = -0.5 * dt * (volatility * volatility * Si * Si / (dS * dS) + r);
            double gamma = 0.25 * dt * (volatility * volatility * Si * Si / (dS * dS) + r * Si / dS);

            if(j !=1){
                a_coeff[j-2] = -alpha; // a1 to a_{n-1}
            }
            else{
                alpha_first = alpha; // Store alpha for the first interior node
            }

            b_coeff[j-1] = 1.0 - beta; // b1 to bn

            if(j !=n){
                c_coeff[j-1] = -gamma; // c1 to c_{n-1}
            }
            else{
                gamma_last = gamma; // Store gamma for the last interior node
            }

            // Right-hand side
            d_coeff[j-1] = alpha * V_prev[j-1] + (1.0 + beta) * V_prev[j] + gamma * V_prev[j+1];
        }

        // Adjust RHS for boundary conditions
        // Lower boundary (S=0)
        if(type == "Call"){
            // V(0, t) = 0
            // d_coeff[0] += alpha_first * V(0, t) = alpha_first * 0 = 0 (no change)
            // No action needed
        }
        else if(type == "Put"){
            // V(0, t) = K * e^{-r*(T-t)}
            double boundary = strike * std::exp(-r * (maturity - step * dt));
            d_coeff[0] += alpha_first * boundary;
        }

        // Upper boundary (S=S_max)
        if(type == "Call"){
            // V(S_max, t) = S_max - K * e^{-r*(T-t)}
            double boundary = S_max - strike * std::exp(-r * (maturity - step * dt));
            d_coeff[n -1] += gamma_last * boundary;
        }
        else if(type == "Put"){
            // V(S_max, t) = 0
            // d_coeff[n -1] += gamma_last * V(S_max, t) = gamma_last * 0 = 0 (no change)
            // No action needed
        }

        // Solve the tridiagonal system
        std::vector<double> solution = Utils::thomasAlgorithm(a_coeff, b_coeff, c_coeff, d_coeff);

        // Check for NaN or Inf in the solution
        for(int j =0; j < solution.size(); ++j){
            if(std::isnan(solution[j]) || std::isinf(solution[j])){
                std::cerr << "Error: Solution contains NaN or Inf at node " << j+1 << "\n";
                exit(EXIT_FAILURE);
            }
        }

        // Update V_new with the solution
        for(int j =1; j <=n; ++j){
            V_new[j] = solution[j-1];
        }

        // Apply boundary conditions
        if(type == "Call"){
            V_new[0] = 0.0;
            V_new[M] = S_max - strike * std::exp(-r * (maturity - step * dt));
        }
        else if(type == "Put"){
            V_new[0] = strike * std::exp(-r * (maturity - step * dt));
            V_new[M] = 0.0;
        }

        // Update V_prev for the next time step
        V_prev = V_new;
    }

    // Interpolate to find the option price at spotPrice
    // Find the closest asset price indices
    if(spotPrice >= S[M]){
        return V_prev[M];
    }
    if(spotPrice <= S[0]){
        return V_prev[0];
    }

    int idx = std::lower_bound(S.begin(), S.end(), spotPrice) - S.begin();
    if(idx == 0){
        return V_prev[0];
    }
    double S1 = S[idx-1];
    double S2 = S[idx];
    double V1 = V_prev[idx-1];
    double V2 = V_prev[idx];

    // Linear interpolation
    double optionPrice = V1 + (V2 - V1) * (spotPrice - S1) / (S2 - S1);

    return optionPrice;
}

// Implement Greek calculations
double EuropeanOption::delta(const Market& market, double dS) const {
    // Clone with perturbed spot price
    EuropeanOption* option_plus = cloneWithPerturbedSpot(spotPrice + dS);
    double price_plus = option_plus->price(market);
    
    EuropeanOption* option_minus = cloneWithPerturbedSpot(spotPrice - dS);
    double price_minus = option_minus->price(market);
    
    // Delta approximation
    double delta = (price_plus - price_minus) / (2 * dS);
    
    // Clean up
    delete option_plus;
    delete option_minus;
    
    return delta;
}

double EuropeanOption::gamma(const Market& market, double dS) const {
    // Clone with perturbed spot prices
    EuropeanOption* option_plus = cloneWithPerturbedSpot(spotPrice + dS);
    double price_plus = option_plus->price(market);
    
    EuropeanOption* option_minus = cloneWithPerturbedSpot(spotPrice - dS);
    double price_minus = option_minus->price(market);
    
    double price_current = this->price(market);
    
    // Gamma approximation
    double gamma = (price_plus - 2 * price_current + price_minus) / (dS * dS);
    
    // Clean up
    delete option_plus;
    delete option_minus;
    
    return gamma;
}

double EuropeanOption::rho(const Market& market, double dr) const {
    // Clone Market with perturbed rates
    Market market_plus = market;
    double currentRate = market_plus.getRate(maturity);
    market_plus.setRate(maturity, currentRate + dr);
    double price_plus = this->price(market_plus);
    
    Market market_minus = market;
    market_minus.setRate(maturity, currentRate - dr);
    double price_minus = this->price(market_minus);
    
    // Rho approximation
    double rho = (price_plus - price_minus) / (2 * dr);
    
    return rho;
}

double EuropeanOption::vega(const Market& market, double dSigma) const {
    // Clone with perturbed volatility
    EuropeanOption* option_plus = cloneWithPerturbedVolatility(volatility + dSigma);
    double price_plus = option_plus->price(market);
    
    EuropeanOption* option_minus = cloneWithPerturbedVolatility(volatility - dSigma);
    double price_minus = option_minus->price(market);
    
    // Vega approximation
    double vega = (price_plus - price_minus) / (2 * dSigma);
    
    // Clean up
    delete option_plus;
    delete option_minus;
    
    return vega;
}

double EuropeanOption::theta(const Market& market, double dT) const {
    if (maturity <= dT) {
        throw std::invalid_argument("Maturity too small to compute Theta.");
    }
    
    // Clone with perturbed maturity
    EuropeanOption* option_new = cloneWithPerturbedMaturity(maturity - dT);
    double price_new = option_new->price(market);
    
    double price_current = this->price(market);
    
    // Theta approximation
    double theta = (price_new - price_current) / dT;
    
    // Clean up
    delete option_new;
    
    return theta;
}
