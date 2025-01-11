# Vanilla Options Pricing Using Finite Difference Methods

## Authors
This project was developed by **Alaa Bouattour** and **Mahdi Ben Ayed** as part of the **M2 Ingénierie et Finance** program.

## Project Overview
The goal of this project is to price vanilla options (European and American) using a numerical solution to the Black–Scholes Partial Differential Equation (PDE). The **Crank–Nicolson finite difference method** is used for the computations, which are implemented in **C++**. 

An **Excel file** is provided as a user interface to control the pricer, allowing users to input parameters, run the model, and visualize results.

## Features

### Inputs
The following inputs are required:
1. **Contract type**: Call or Put.
2. **Exercise type**: European or American.
3. **Maturity**: Time to expiration (in years).
4. **Strike price**: The exercise price of the option.
5. **Calculation date**: Defaults to the current date if not provided.
6. **Time discretization**: Number of time steps for the finite difference grid.
7. **Spot discretization**: Number of spot steps for the finite difference grid.
8. **Current spot price (S₀)**: Price of the underlying asset.
9. **Risk-free rate**: A piecewise linear function defined as pairs \((t_k, r(t_k))\).
10. **Volatility (σ)**: Annualized volatility of the underlying asset.

### Outputs
The program produces the following outputs:
1. **Option price** (\(P\)).
2. **Greek sensitivities**:
   - Delta (\(Δ\)): Sensitivity of \(P\) to changes in \(S\).
   - Gamma (\(Γ\)): Sensitivity of \(Δ\) to changes in \(S\).
   - Theta (\(Θ\)): Sensitivity of \(P\) to changes in \(t\).
   - Rho (\(ρ\)): Sensitivity of \(P\) to changes in \(r\).
   - Vega (\(v\)): Sensitivity of \(P\) to changes in \(σ\).
3. **Graphs**:
   - Option price (\(P(S, T₀)\)) as a function of the underlying price.
   - Delta (\(Δ(S, T₀)\)) as a function of the underlying price.
   - Exercise boundary for American options.

### Validation
Results are validated using:
- The **Black–Scholes formula** for European options.
- Comparisons between American and European options in scenarios where their behavior should match.

## File Structure

### Source Code (C++)
- Implements the Crank–Nicolson finite difference scheme.
- Outputs the results as CSV files.

### Excel Interface
The Excel file (`.xlsm`) serves as the **user interface**. It:
1. Collects input parameters from the user.
2. Executes the C++ pricer through macros.
3. Loads and displays results (option prices, Greeks, and graphs).

### Output Files
- `option_greeks.csv`: Contains the Greek sensitivities.
- `option_price_vs_S.csv`: Contains the option prices as a function of the underlying price.
- `delta_vs_S.csv`: Contains delta values as a function of the underlying price.

## Usage Instructions

### Prerequisites
1. Install a C++ compiler (e.g., GCC).
2. Ensure Excel supports macros and enable them in the settings.

### Steps
1. Open the **Excel file**.
2. Fill in the required inputs in the `Parameters` sheet.
3. Run the pricer by clicking the appropriate macro button.
4. Results will be displayed in dedicated sheets or loaded from the output CSV files.

### Compilation and Execution
The pricer can be compiled manually using:
```bash
g++ -std=c++17 main.cpp EuropeanOption.cpp AmericanOption.cpp Market.cpp Utils.cpp -o option_pricer  
./option_pricer
