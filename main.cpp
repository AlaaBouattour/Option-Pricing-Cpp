// main.cpp
#include "EuropeanOption.hpp"
#include "AmericanOption.hpp"
#include "Market.hpp"
#include "Option.hpp"
#include "Utils.hpp"
#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>

int main() {
    try {
        // File names
        std::string ratesFile = "inputs_rates.csv";      // Interest rates
        std::string optionsFile = "inputs_model.csv";  // Option parameters

        // Create and load Market instance
        Market market;
        market.loadRates(ratesFile); // Ensure this file exists and is correctly formatted
        std::cout << "[Market] Rates loaded from " << ratesFile << ".\n";
        market.displayRates();

        // Open the options file
        std::ifstream infile(optionsFile);
        if (!infile.is_open()) {
            throw std::runtime_error("Unable to open options file: " + optionsFile);
        }

        std::string line;
        // Read the header line
        if (!std::getline(infile, line)) {
            throw std::runtime_error("Options file is empty.");
        }

        // Create a map to hold parameter-value pairs
        std::map<std::string, std::string> paramMap;

        // Read each line and populate the map
        while (std::getline(infile, line)) {
            if (line.empty()) continue; // Skip empty lines

            std::stringstream ss(line);
            std::string param, value;

            // Read until the first comma for parameter
            if (!std::getline(ss, param, ',')) continue;
            // Read the rest for value (in case value contains commas)
            if (!std::getline(ss, value)) continue;

            // Trim whitespace from param and value
            param.erase(0, param.find_first_not_of(" \t\r\n")); // Left trim
            param.erase(param.find_last_not_of(" \t\r\n") + 1); // Right trim
            value.erase(0, value.find_first_not_of(" \t\r\n")); // Left trim
            value.erase(value.find_last_not_of(" \t\r\n") + 1); // Right trim

            paramMap[param] = value;
        }

        infile.close();

        // Extract parameters from the map
        std::string type = paramMap["Type de contrat"];           // "Call" or "Put"
        std::string exerciseType = paramMap["Type d'exercice"];  // "européen" or "américain"
        double maturity = std::stod(paramMap["Maturité (en années)"]);
        double strike = std::stod(paramMap["Prix d’exercice (strike)"]);
        double spotPrice = std::stod(paramMap["Prix actuel (S0)"]);
        double volatility = std::stod(paramMap["Volatilité (sigma)"]);
        int N = std::stoi(paramMap["Discrétisation (temps)"]);   // Time steps
        int M = std::stoi(paramMap["Discrétisation (spot)"]);    // Asset price steps
        std::string calculationDate = paramMap["Date de calcul"];

        // Create Option instance using smart pointers
        std::unique_ptr<Option> myOption;

        if (exerciseType == "européen") {
            myOption = std::make_unique<EuropeanOption>(
                type, exerciseType, maturity, strike,
                calculationDate, spotPrice, volatility, M, N
            );
        }
        else if (exerciseType == "américain") {
            myOption = std::make_unique<AmericanOption>(
                type, exerciseType, maturity, strike,
                calculationDate, spotPrice, volatility, M, N
            );
        }
        else {
            std::cerr << "Invalid exercise type specified.\n";
            return 1;
        }

        // Display option details
        myOption->display();

        // Calculate Option Price
        double optionPrice = myOption->price(market);
        std::cout << "Option Price: " << optionPrice << "\n";

        // Calculate Greeks
        // Downcast to access derived class methods
        if (exerciseType == "européen") {
            EuropeanOption* euroOpt = dynamic_cast<EuropeanOption*>(myOption.get());
            if (euroOpt) {
                double dS = 1.0;         // $1
                double delta = euroOpt->delta(market, dS);
                double gamma = euroOpt->gamma(market, dS);
                double dr = 0.0001;      // 0.01%
                double rho = euroOpt->rho(market, dr);
                double dSigma = 0.0001;  // 0.01%
                double vega = euroOpt->vega(market, dSigma);
                double dT = 1e-4;        // ~0.0365 days
                double theta = euroOpt->theta(market, dT);

                std::cout << "Delta: " << delta << "\n";
                std::cout << "Gamma: " << gamma << "\n";
                std::cout << "Rho: " << rho << "\n";
                std::cout << "Vega: " << vega << "\n";
                std::cout << "Theta: " << theta << "\n";
            }
            else {
                std::cerr << "Error: Failed to cast to EuropeanOption.\n";
            }
        }
        else if (exerciseType == "américain") {
            AmericanOption* amOpt = dynamic_cast<AmericanOption*>(myOption.get());
            if (amOpt) {
                double dS = 1.0;         // $1
                double delta = amOpt->delta(market, dS);
                double gamma = amOpt->gamma(market, dS);
                double dr = 0.0001;      // 0.01%
                double rho = amOpt->rho(market, dr);
                double dSigma = 0.0001;  // 0.01%
                double vega = amOpt->vega(market, dSigma);
                double dT = 1e-4;        // ~0.0365 days
                double theta = amOpt->theta(market, dT);

                std::cout << "Delta: " << delta << "\n";
                std::cout << "Gamma: " << gamma << "\n";
                std::cout << "Rho: " << rho << "\n";
                std::cout << "Vega: " << vega << "\n";
                std::cout << "Theta: " << theta << "\n";
            }
            else {
                std::cerr << "Error: Failed to cast to AmericanOption.\n";
            }
        }

    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
