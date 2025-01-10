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
            param.erase(0, param.find_first_not_of(" \t\r\n")); // Suppression gauche
            param.erase(param.find_last_not_of(" \t\r\n") + 1); // Suppression droite
            value.erase(0, value.find_first_not_of(" \t\r\n")); // Suppression gauche
            value.erase(value.find_last_not_of(" \t\r\n") + 1); // Suppression droite

            paramMap[param] = value;
        }

        infile.close();

        // Extraction des paramètres de la map
        std::string type = paramMap["Type de contrat"];           // "Call" ou "Put"
        std::string exerciseType = paramMap["Type d exercice"];  // "européen" ou "américain"
        double maturity = std::stod(paramMap["Maturite (en annees)"]);
        double strike = std::stod(paramMap["Prix d exercice (strike)"]);
        double spotPrice = std::stod(paramMap["Prix actuel (S0)"]);
        double volatility = std::stod(paramMap["Volatilite (sigma)"]);
        int N = std::stoi(paramMap["Discretisation (temps)"]);   // Pas de temps
        int M = std::stoi(paramMap["Discretisation (spot)"]);    // Pas de prix de l'actif
        std::string calculationDate = paramMap["Date de calcul"];

        // Création de l'instance Option en utilisant des pointeurs intelligents
        std::unique_ptr<Option> myOption;

        if (exerciseType == "europeen") {
            myOption = std::make_unique<EuropeanOption>(
                type, exerciseType, maturity, strike,
                calculationDate, spotPrice, volatility, M, N
            );
        }
        else if (exerciseType == "americain") {
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

        // Déclaration des variables pour les Grecques
        double delta = 0.0;
        double gamma = 0.0;
        double rho = 0.0;
        double vega = 0.0;
        double theta = 0.0;

        // Calcul des Grecques
        // Downcast pour accéder aux méthodes des classes dérivées
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
        else if (exerciseType == "americain") {
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

        // -------------------------------
        // EXPORTATION DU PRIX DE L'OPTION ET DES GRECQUES
        // -------------------------------

        // nouveau fichier CSV pour stocker le prix de l'option et les Grecques
        std::ofstream greeksFile("option_greeks.csv");
        if (!greeksFile.is_open()) {
            throw std::runtime_error("Impossible de créer le fichier de sortie pour les grecques.");
        }

        // l'en-tête
        greeksFile << "Option_Price,Delta,Gamma,Rho,Vega,Theta\n";

        greeksFile << optionPrice << "," << delta << "," << gamma << "," << rho << "," << vega << "," << theta << "\n";

        greeksFile.close();

        std::cout << "[Données] Le prix de l'option et les Grecques ont été écrits dans 'option_greeks.csv'.\n";

        // -------------------------------
        // EXPORTATION DU PRIX DE L'OPTION EN FONCTION DE S
        // -------------------------------

        // plage de valeurs S pour le tracé
        double S_min = spotPrice * 0.5; // 50% du prix initial
        double S_max = spotPrice * 1.5; // 150% du prix initial
        int numPoints = 100;             // Nombre de points dans le tracé

        double delta_S = (S_max - S_min) / (numPoints - 1);
        double delta_S_finite = 1.0;     // Pas pour le calcul des Grecques (∆)

        // vecteurs pour stocker S et P(S,T)
        std::vector<double> S_values;
        std::vector<double> P_values;

        // calcul de P(S,T)
        for(int i = 0; i < numPoints; ++i){
            double current_S = S_min + i * delta_S;

            // instance temporaire de l'Option avec current_S
            std::unique_ptr<Option> tempOption;

            if (exerciseType == "europeen") {
                tempOption = std::make_unique<EuropeanOption>(
                    type, exerciseType, maturity, strike,
                    calculationDate, current_S, volatility, M, N
                );
            }
            else if (exerciseType == "americain") {
                tempOption = std::make_unique<AmericanOption>(
                    type, exerciseType, maturity, strike,
                    calculationDate, current_S, volatility, M, N
                );
            }

            // prix de l'option pour current_S
            double price = tempOption->price(market);

            // Stocker les valeurs
            S_values.push_back(current_S);
            P_values.push_back(price);
        }

        // S et P(S,T) --> CSV
        std::ofstream outfile("option_price_vs_S.csv");
        if (!outfile.is_open()) {
            throw std::runtime_error("Impossible de créer le fichier de sortie pour le tracé.");
        }

        // l'en-tête
        outfile << "S,P(S,T)\n";

        for(int i =0; i < numPoints; ++i){
            outfile << S_values[i] << "," << P_values[i] << "\n";
        }

        outfile.close();

        std::cout << "\n[Données] Les prix de l'option en fonction de S ont été écrits dans 'option_price_vs_S.csv'.\n";

        // -------------------------------
        // CALCUL ET EXPORTATION DE DELTA(S,T0)
        // -------------------------------

        // vecteur pour stocker Delta(S,T0)
        std::vector<double> Delta_values;

        // Delta(S,T0) avec les différences finies
        for(int i = 0; i < numPoints; ++i){
            if(i == 0){
                // Différence avant
                double P_plus = P_values[i+1];
                double P_current = P_values[i];
                double delta_S_local = S_values[i+1] - S_values[i];
                double Delta = (P_plus - P_current) / delta_S_local;
                Delta_values.push_back(Delta);
            }
            else if(i == numPoints -1){
                // Différence arrière
                double P_current = P_values[i];
                double P_minus = P_values[i-1];
                double delta_S_local = S_values[i] - S_values[i-1];
                double Delta = (P_current - P_minus) / delta_S_local;
                Delta_values.push_back(Delta);
            }
            else{
                // Différence centrale
                double P_plus = P_values[i+1];
                double P_minus = P_values[i-1];
                double Delta = (P_plus - P_minus) / (2 * delta_S_finite);
                Delta_values.push_back(Delta);
            }
        }

        // Delta(S,T0) --> CSV
        std::ofstream deltaFile("delta_vs_S.csv");
        if (!deltaFile.is_open()) {
            throw std::runtime_error("Impossible de créer le fichier de sortie pour Delta.");
        }

        // l'en-tête
        deltaFile << "S,Delta(S,T0)\n";

        for(int i =0; i < numPoints; ++i){
            deltaFile << S_values[i] << "," << Delta_values[i] << "\n";
        }

        deltaFile.close();

        std::cout << "[Données] Delta(S,T0) en fonction de S ont été écrits dans 'delta_vs_S.csv'.\n";

    }
    catch (const std::exception& e) {
        std::cerr << "Exception : " << e.what() << "\n";
        return 1;
    }

    return 0;
}
