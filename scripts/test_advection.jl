using ExoClouds.Physics: dmp_dt
using Unitful
using Unitful: km, m, s
using DifferentialEquations

N = 30
el = water()
atm = Atmosphere(71492km, 24.79m/s^2, (1:N) * km, fill(2.2u, N), 10 .^(-5:8/(N-1):3) * bar)
