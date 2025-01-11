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
#include <vector>

int main() {
    try {
        // Noms des fichiers
        std::string ratesFile = "inputs_rates.csv";      // Taux d'intérêt
        std::string optionsFile = "inputs_model.csv";    // Paramètres de l'option

        // Création et chargement de l'instance Market
        Market market;
        market.loadRates(ratesFile); 
        std::cout << "[Market] Taux chargés depuis " << ratesFile << ".\n";
        market.displayRates();

        // Ouverture du fichier des options
        std::ifstream infile(optionsFile);
        if (!infile.is_open()) {
            throw std::runtime_error("Impossible d'ouvrir le fichier des options : " + optionsFile);
        }

        std::string line;
        // Lecture de la ligne d'en-tête
        if (!std::getline(infile, line)) {
            throw std::runtime_error("Le fichier des options est vide.");
        }

        // Création d'une map pour stocker les paires paramètre-valeur
        std::map<std::string, std::string> paramMap;

        // Lecture de chaque ligne et remplissage de la map
        while (std::getline(infile, line)) {
            if (line.empty()) continue; // Ignorer les lignes vides

            std::stringstream ss(line);
            std::string param, value;

            // Lecture jusqu'à la première virgule pour le paramètre
            if (!std::getline(ss, param, ',')) continue;
            // Lecture du reste pour la valeur (au cas où la valeur contiendrait des virgules)
            if (!std::getline(ss, value)) continue;

            // Suppression des espaces blancs au début et à la fin
            param.erase(0, param.find_first_not_of(" \t\r\n")); 
            param.erase(param.find_last_not_of(" \t\r\n") + 1); 
            value.erase(0, value.find_first_not_of(" \t\r\n")); 
            value.erase(value.find_last_not_of(" \t\r\n") + 1); 

            paramMap[param] = value;
        }

        infile.close();

        // Extraction des paramètres de la map
        std::string type = paramMap["Type de contrat"];            // "Call" ou "Put"
        std::string exerciseType = paramMap["Type d exercice"];    // "europeen" ou "american"
        double maturity = std::stod(paramMap["Maturite (en annees)"]);
        double strike = std::stod(paramMap["Prix d exercice (strike)"]);
        double spotPrice = std::stod(paramMap["Prix actuel (S0)"]);
        double volatility = std::stod(paramMap["Volatilite (sigma)"]);
        int N = std::stoi(paramMap["Discretisation (temps)"]);     // Pas de temps
        int M = std::stoi(paramMap["Discretisation (spot)"]);      // Pas de prix de l'actif
        std::string calculationDate = paramMap["Date de calcul"];

        // Création de l'instance Option en utilisant des pointeurs intelligents
        std::unique_ptr<Option> myOption;

        if (exerciseType == "europeen") {
            myOption = std::make_unique<EuropeanOption>(
                type, exerciseType, maturity, strike,
                calculationDate, spotPrice, volatility, M, N
            );
        }
        else if (exerciseType == "american") {
            myOption = std::make_unique<AmericanOption>(
                type, exerciseType, maturity, strike,
                calculationDate, spotPrice, volatility, M, N
            );
        }
        else {
            std::cerr << "Type d'exercice invalide spécifié.\n";
            return 1;
        }

        // Affichage des détails de l'option
        myOption->display();

        // Calcul du prix de l'option au prix initial de l'actif
        double optionPrice = myOption->price(market);
        std::cout << "Option Price: " << optionPrice << "\n";

        // Déclaration des variables pour les Grecques globales
        double delta = 0.0;
        double gamma = 0.0;
        double rho = 0.0;
        double vega = 0.0;
        double theta = 0.0;

        // Calcul des Grecques (globales) via downcast
        if (exerciseType == "europeen") {
            EuropeanOption* euroOpt = dynamic_cast<EuropeanOption*>(myOption.get());
            if (euroOpt) {
                double dS = 1.0;         // $1
                delta = euroOpt->delta(market, dS);
                gamma = euroOpt->gamma(market, dS);
                double dr = 0.0001;      // 0.01%
                rho = euroOpt->rho(market, dr);
                double dSigma = 0.0001;  // 0.01%
                vega = euroOpt->vega(market, dSigma);
                double dT = 1e-4;        // ~0.0365 jours
                theta = euroOpt->theta(market, dT);

                std::cout << "Delta: " << delta << "\n";
                std::cout << "Gamma: " << gamma << "\n";
                std::cout << "Rho: " << rho << "\n";
                std::cout << "Vega: " << vega << "\n";
                std::cout << "Theta: " << theta << "\n";
            }
            else {
                std::cerr << "Erreur : Échec du cast en EuropeanOption.\n";
            }
        }
        else if (exerciseType == "american") {
            AmericanOption* amOpt = dynamic_cast<AmericanOption*>(myOption.get());
            if (amOpt) {
                double dS = 1.0;         // $1
                delta = amOpt->delta(market, dS);
                gamma = amOpt->gamma(market, dS);
                double dr = 0.0001;      // 0.01%
                rho = amOpt->rho(market, dr);
                double dSigma = 0.0001;  // 0.01%
                vega = amOpt->vega(market, dSigma);
                double dT = 1e-4;        // ~0.0365 jours
                theta = amOpt->theta(market, dT);

                std::cout << "Delta: " << delta << "\n";
                std::cout << "Gamma: " << gamma << "\n";
                std::cout << "Rho: " << rho << "\n";
                std::cout << "Vega: " << vega << "\n";
                std::cout << "Theta: " << theta << "\n";
            }
            else {
                std::cerr << "Erreur : Échec du cast en AmericanOption.\n";
            }
        }

        // Exportation globale du prix et des Grecques
        std::ofstream greeksFile("option_greeks.csv");
        if (!greeksFile.is_open()) {
            throw std::runtime_error("Impossible de créer le fichier de sortie pour les grecques.");
        }

        greeksFile << "Option_Price,Delta,Gamma,Rho,Vega,Theta\n";
        greeksFile << optionPrice << "," << delta << "," << gamma << "," << rho << "," << vega << "," << theta << "\n";
        greeksFile.close();

        std::cout << "[Données] Le prix de l'option et les Grecques ont été écrits dans 'option_greeks.csv'.\n";

        // -------------------------------
        // EXPORTATION DU PRIX DE L'OPTION EN FONCTION DE S
        // -------------------------------
        double S_min = strike * 0.3; 
        double S_max = strike * 2.0; 
        int numPoints = 100; 

        double delta_S = (S_max - S_min) / (numPoints - 1);

        // Vecteurs pour stocker S et P(S,T)
        std::vector<double> S_values;
        std::vector<double> P_values;

        // Calcul de P(S,T) par instanciation temporaire
        for(int i = 0; i < numPoints; ++i){
            double current_S = S_min + i * delta_S;

            std::unique_ptr<Option> tempOption;
            if (exerciseType == "europeen") {
                tempOption = std::make_unique<EuropeanOption>(
                    type, exerciseType, maturity, strike,
                    calculationDate, current_S, volatility, M, N
                );
            }
            else {
                tempOption = std::make_unique<AmericanOption>(
                    type, exerciseType, maturity, strike,
                    calculationDate, current_S, volatility, M, N
                );
            }

            double price = tempOption->price(market);

            S_values.push_back(current_S);
            P_values.push_back(price);
        }

        // Écriture de P(S,T) dans un CSV
        std::ofstream outfile("option_price_vs_S.csv");
        if (!outfile.is_open()) {
            throw std::runtime_error("Impossible de créer le fichier de sortie pour le tracé.");
        }
        outfile << "S,P(S,T)\n";
        for(int i =0; i < numPoints; ++i){
            outfile << S_values[i] << "," << P_values[i] << "\n";
        }
        outfile.close();

        std::cout << "\n[Données] Les prix de l'option en fonction de S ont été écrits dans 'option_price_vs_S.csv'.\n";

        // -------------------------------
        // CALCUL ET EXPORTATION DE DELTA(S) VIA LA MÉTHODE delta(...)
        // -------------------------------
        double dS_for_delta = 1.0;  // Pas à utiliser pour le calcul de Delta
        std::vector<double> Delta_values_method2;

        for(int i = 0; i < numPoints; ++i){
            double current_S = S_min + i * delta_S;

            // Création d'une option temporaire avec spot = current_S
            std::unique_ptr<Option> tempOption;
            if (exerciseType == "europeen") {
                tempOption = std::make_unique<EuropeanOption>(
                    type, exerciseType, maturity, strike,
                    calculationDate, current_S, volatility, M, N
                );
            }
            else {
                tempOption = std::make_unique<AmericanOption>(
                    type, exerciseType, maturity, strike,
                    calculationDate, current_S, volatility, M, N
                );
            }

            // Calcul de Delta via la méthode delta(...)
            double delta_i = 0.0;
            if (exerciseType == "europeen") {
                EuropeanOption* euroOpt = dynamic_cast<EuropeanOption*>(tempOption.get());
                if(euroOpt) {
                    delta_i = euroOpt->delta(market, dS_for_delta);
                }
            }
            else {
                AmericanOption* amOpt = dynamic_cast<AmericanOption*>(tempOption.get());
                if(amOpt) {
                    delta_i = amOpt->delta(market, dS_for_delta);
                }
            }

            Delta_values_method2.push_back(delta_i);
        }

        // Exporter Delta(S) dans un CSV
        std::ofstream deltaFileMethod2("delta_vs_S.csv");
        if (!deltaFileMethod2.is_open()) {
            throw std::runtime_error("Impossible de créer 'delta_vs_S.csv'.");
        }
        deltaFileMethod2 << "S,Delta(S)\n";
        for(int i =0; i < numPoints; ++i){
            double current_S = S_min + i * delta_S;
            deltaFileMethod2 << current_S << "," << Delta_values_method2[i] << "\n";
        }
        deltaFileMethod2.close();

        std::cout << "[Données] Delta(S) via la méthode delta(...) a été écrit dans 'delta_vs_S.csv'.\n";

    }
    catch (const std::exception& e) {
        std::cerr << "Exception : " << e.what() << "\n";
        return 1;
    }

    return 0;
}
