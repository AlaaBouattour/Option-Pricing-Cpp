// AmericanOption.cpp
#include "AmericanOption.hpp"
#include "Market.hpp"    // Include Market.hpp to access Market class methods
#include "Utils.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>

// Implement clone methods
AmericanOption* AmericanOption::cloneWithPerturbedSpot(double perturbedSpot) const {
    return new AmericanOption(type, exerciseType, maturity, strike,
                              calculationDate, perturbedSpot, volatility, M, N);
}

AmericanOption* AmericanOption::cloneWithPerturbedVolatility(double perturbedVolatility) const {
    return new AmericanOption(type, exerciseType, maturity, strike,
                              calculationDate, spotPrice, perturbedVolatility, M, N);
}

AmericanOption* AmericanOption::cloneWithPerturbedStrike(double perturbedStrike) const {
    return new AmericanOption(type, exerciseType, maturity, perturbedStrike,
                              calculationDate, spotPrice, volatility, M, N);
}

AmericanOption* AmericanOption::cloneWithPerturbedMaturity(double perturbedMaturity) const {
    return new AmericanOption(type, exerciseType, perturbedMaturity, strike,
                              calculationDate, spotPrice, volatility, M, N);
}

// Implement display method
void AmericanOption::display() const {
    std::cout << "American Option Details:\n";
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

// Method to record exercise boundary
void AmericanOption::recordExerciseBoundary(double t, double S_star) const {
    exerciseBoundary_t.push_back(t);
    exerciseBoundary_S_star.push_back(S_star);
}

// Implement price method using Finite Difference Method with Early Exercise
double AmericanOption::price(const Market& market) const {
    // Retrieve the risk-free rate at current time
    double r = market.getRate(maturity);

    // Asset price range
    double S_max = spotPrice * 3.0; // Adjust as needed based on volatility and maturity
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
        double t = step * dt; // Current time

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
            double boundary = strike * std::exp(-r * (maturity - t));
            d_coeff[0] += alpha_first * boundary;
        }

        // Upper boundary (S=S_max)
        if(type == "Call"){
            // V(S_max, t) = S_max - K * e^{-r*(T-t)}
            double boundary = S_max - strike * std::exp(-r * (maturity - t));
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
            V_new[M] = S_max - strike * std::exp(-r * (maturity - t));
        }
        else if(type == "Put"){
            V_new[0] = strike * std::exp(-r * (maturity - t));
            V_new[M] = 0.0;
        }

        // Early exercise condition: compare with intrinsic value
        for(int j =1; j <=n; ++j){
            double intrinsic;
            if(type == "Call"){
                intrinsic = std::max(S[j] - strike, 0.0);
            }
            else if(type == "Put"){
                intrinsic = std::max(strike - S[j], 0.0);
            }
            else{
                throw std::runtime_error("Invalid option type.");
            }
            V_new[j] = std::max(V_new[j], intrinsic);
        }

        // Identify the exercise boundary S_star(t)
        double S_star = 0.0;
        for(int j =1; j <=n; ++j){
            double intrinsic;
            if(type == "Call"){
                intrinsic = std::max(S[j] - strike, 0.0);
            }
            else if(type == "Put"){
                intrinsic = std::max(strike - S[j], 0.0);
            }

            if(V_new[j] == intrinsic){
                S_star = S[j];
                break; // Found the critical point
            }
        }

        // Record the exercise boundary
        recordExerciseBoundary(t, S_star);

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

    // Export the exercise boundary to CSV
    std::ofstream boundaryFile("exercise_boundary.csv");
    if (!boundaryFile.is_open()) {
        throw std::runtime_error("Unable to create 'exercise_boundary.csv'.");
    }

    // Write header
    boundaryFile << "t,S_star(t)\n";

    // Write data
    for(int i =0; i < exerciseBoundary_t.size(); ++i){
        boundaryFile << exerciseBoundary_t[i] << "," << exerciseBoundary_S_star[i] << "\n";
    }

    boundaryFile.close();

    std::cout << "[Data] Exercise boundary has been written to 'exercise_boundary.csv'.\n";

    return optionPrice;
}

// Implement Greek calculations
double AmericanOption::delta(const Market& market, double dS) const {
    // Clone with perturbed spot price
    AmericanOption* option_plus = cloneWithPerturbedSpot(spotPrice + dS);
    double price_plus = option_plus->price(market);
    
    AmericanOption* option_minus = cloneWithPerturbedSpot(spotPrice - dS);
    double price_minus = option_minus->price(market);
    
    // Delta approximation
    double delta = (price_plus - price_minus) / (2 * dS);
    
    // Clean up
    delete option_plus;
    delete option_minus;
    
    return delta;
}

double AmericanOption::gamma(const Market& market, double dS) const {
    // Clone with perturbed spot prices
    AmericanOption* option_plus = cloneWithPerturbedSpot(spotPrice + dS);
    double price_plus = option_plus->price(market);
    
    AmericanOption* option_minus = cloneWithPerturbedSpot(spotPrice - dS);
    double price_minus = option_minus->price(market);
    
    double price_current = this->price(market);
    
    // Gamma approximation
    double gamma = (price_plus - 2 * price_current + price_minus) / (dS * dS);
    
    // Clean up
    delete option_plus;
    delete option_minus;
    
    return gamma;
}

double AmericanOption::rho(const Market& market, double dr) const {
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

double AmericanOption::vega(const Market& market, double dSigma) const {
    // Clone with perturbed volatility
    AmericanOption* option_plus = cloneWithPerturbedVolatility(volatility + dSigma);
    double price_plus = option_plus->price(market);
    
    AmericanOption* option_minus = cloneWithPerturbedVolatility(volatility - dSigma);
    double price_minus = option_minus->price(market);
    
    // Vega approximation
    double vega = (price_plus - price_minus) / (2 * dSigma);
    
    // Clean up
    delete option_plus;
    delete option_minus;
    
    return vega;
}

double AmericanOption::theta(const Market& market, double dT) const {
    if (maturity <= dT) {
        throw std::invalid_argument("Maturity too small to compute Theta.");
    }
    
    // Clone with perturbed maturity
    AmericanOption* option_new = cloneWithPerturbedMaturity(maturity - dT);
    double price_new = option_new->price(market);
    
    double price_current = this->price(market);
    
    // Theta approximation
    double theta = (price_new - price_current) / dT;
    
    // Clean up
    delete option_new;
    
    return theta;
}
