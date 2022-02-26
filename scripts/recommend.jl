using ExoClouds
using Unitful: bar, K, u

pressure = 10 .^(-5:8/29:3) * bar
temperature = fill(1300.0, 30) * K
metallicity = 1
mean_molecular_weight = 2.2u

recommended = recommend_gas(pressure, temperature, metallicity, mean_molecular_weight)
