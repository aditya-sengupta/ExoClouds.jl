# just put in all the equations from the appendix of Gao Ackerman Marley (2018)

using Unitful
using Unitful: k, R

function critical_radius(M::Mass, σₛ::Energy, ρₚ::Density, S::Real, T::Temperature)
    return 2 * M * σₛ / (ρₚ * R * T * log(S))
end

function diffusion_rate(e::Element, p::Pressure, T::Temperature) # p for partial pressure
    return p * mixing_ratio(m) / sqrt(2π * molecular_weight(m) * k * T)
end

function zeldovich(F, T, gₘ)
    return sqrt(F / (3π * k * T * gₘ^2))
end

"""
Homogeneous nucleation.

m - the element nucleating

"""
function nucleation(e::Element, σₛ::Energy, ρₚ::Density, S::Real, p::Pressure, T::Temperature, n::Integer, gₘ::Integer)
    a = critical_radius(molecular_weight(m), σₛ, ρₚ, S, T)
    F = (4π/3) * σₛ * a^2
    Φ = diffusion_rate(m, p, T)
    Z = zeldovich(F, T, gₘ)
    return 4π * a^2 * Φ * Z * n * exp(-F / (k * T))
end

function nucleation(agent::Element, surface::Element, ..., p::Pressure, T::Temperature, gₘ::Integer) 
    Φ = diffusion_rate(agent, p, T)
    c_surf = (Φ / ν) * exp(F_des / (k * T))
    Z = zeldovich(F, T, gₘ)
    x = r / a
    ϕ = sqrt(1 - 2 * μ * x + x^2)
    f₀ = (x - μ) / ϕ
    f = (1/2) * (1 + ((1 - μ * x) / ϕ)^3 + x^3 * (2 - 3 * f₀ + f₀^3) + 3 * μ * x^2 (f₀ - 1))
    return 4π^2 * r_CN^2 * Φ * c_surf * Z * exp(-F*f / (k * T))
end