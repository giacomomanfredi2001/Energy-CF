# Energy-CF

## Overview of the project:

This project explores the modeling of French electricity futures using two distinct stochastic approaches to capture key market dynamics:

OU-IG Model: An Arithmetic Ornstein-Uhlenbeck process with Inverse Gaussian jumps, incorporating mean reversion, seasonality, and price jumps.

OU-Normal Model: A simplified Gaussian-based Ornstein-Uhlenbeck model, designed for parsimony while maintaining essential price dynamics.

## Key Components:

Spot Price Dynamics – Integrating seasonality and stochastic processes.

Forward and Swap Pricing – Deriving analytical solutions under a risk-neutral measure.

Model Calibration – Using optimization techniques to fit market data.

Monte Carlo Simulation – Simulating option prices and evaluating model robustness.

The dataset consists of French power futures with various delivery periods, allowing for a comparative analysis of model complexity and its implications for energy markets.

## Files and Resources

MATLAB Code – Scripts for model calibration, pricing, and simulations.

Project Report (Energy_Finance_Report.pdf) – A detailed discussion of methodologies, results, and model comparisons.

Market Data – Prices of French electricity futures and relevant inputs for calibration.
