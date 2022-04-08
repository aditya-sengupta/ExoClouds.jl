# just put in all the equations from the appendix of Gao Ackerman Marley (2018)

using Unitful
using Unitful: k, R

function critical_radius(M::Mass, σₛ, ρₚ, S, T::Temperature)
    return 2 * M * σₛ / (ρₚ * R * T * log(S))
end

function diffusion_rate(p::Pressure, m::Mass, T::Temperature) # p for partial pressure
    return p / sqrt(2π * m * k * T)
end

function zeldovich(F, T, gₘ)
    return sqrt(F / (3π * k * T * gₘ^2))
end

function homogeneous_nuc_rate(M::Mass, σₛ::Energy, ρₚ::Density, S::Real, T::Temperature, n, Z,)
    a = critical_radius(M, σₛ, ρₚ, S, T)
    F = (4π/3) * σₛ * a^2
    Φ = diffusion_rate(p, m, T) # what is m? need to loop this by particle?
    Z = zeldovich(F, T, gₘ) # what is gₘ
    return 4π * a^2 * Φ * Z * n * exp(-F / (k * T))
end

function heterogeneous_nuc_rate() # to figure this out
    Φ = diffusion_rate(p, m, T) # what is m? need to loop this by particle?
    c_surf = (Φ / ν) * exp(F_des / (k * T))
    Z = zeldovich(F, T, gₘ) # what is gₘ
    x = r / a
    ϕ = sqrt(1 - 2 * μ * x + x^2)
    f₀ = (x - μ) / ϕ
    f = (1/2) * (1 + ((1 - μ * x) / ϕ)^3 + x^3 * (2 - 3 * f₀ + f₀^3) + 3 * μ * x^2 (f₀ - 1))
    return 4π^2 * r_CN^2 * Φ * c_surf * Z * exp(-F*f / (k * T))
end