# just put in all the equations from the appendix of Gao Ackerman Marley (2018)

using Unitful
using Unitful: k, R
using Unitful: erg, Hz

function critical_radius(M::Mass, σₛ::Energy, ρₚ::Density, S::Float64, T::Temperature)
    return 2 * M * σₛ / (ρₚ * R * T * log(S))
end

"""
p is the partial pressure
"""
function diffusion_rate(e::Element, p::Pressure, T::Temperature)
    return p / sqrt(2π * molecular_weight(e) * k * T)
end

function zeldovich(F, T, gₘ)
    return sqrt(F / (3π * k * T * gₘ^2))
end

"""
Homogeneous nucleation (homnucgen.F90)

e - the element nucleating and condensing
n - number density of condensate 
gₘ - # molecules in particles with radius = the critical radius
"""
function hom_nucleation(e::Element, gas_conc::Density, T::Temperature; is_ice=false)
    S = supersaturation(e, gas_conc, T; is_ice=is_ice)
    ρ, m = density(e), molecular_weight(e)
    a = critical_radius(m, σₛ, ρ, S, T)
    σₛ = surface_tension(e, T)
    F = (4π/3) * σₛ * a^2
    Φ = diffusion_rate(e, gas_conc * R * T, T)
    Z = zeldovich(F, T, 4π/3 * ρ * a^3 / m)
    n = gas_conc / molecular_weight(e)
    return 4π * a^2 * Φ * Z * n * exp(-F / (k * T))
end

"""
The heterogeneous nucleation rate is independent of the element being nucleated onto, other than via the radius r = r(ibin,igroup). So we specify the radius as a length alongside the particle it's associated with.

I'm completely ignoring the CARMA implementation here and just going off Pruppacher and Klett (except for where they set mu = 0.95)
"""
function het_nucleation(
    e::Element, r::Length, gas_conc::Density, T::Temperature; 
    μ::Float64=0.95, F_des::FloatEner=2.9e-13*erg, ν::FloatFreq=1e13*Hz, is_ice=false
)
    S = supersaturation(e, gas_conc, T; is_ice=is_ice)
    ρ, M = density(e), molecular_weight(e)
    σₛ = surface_tension(e, T)
    a = critical_radius(M, σₛ, ρ, S, T)
    Φ = diffusion_rate(e, gas_conc * R * T, T)
    F = (4π/3) * σₛ * a^2
    x = r / a
    ϕ = sqrt(1 - 2 * μ * x + x^2)
    f₀ = (x - μ) / ϕ
    f = (1/2) * (1 + ((1 - μ * x) / ϕ)^3 + x^3 * (2 - 3 * f₀ + f₀^3) + 3 * μ * x^2 * (f₀ - 1))
    exp_argument = (F_des - F * f) / k * T
    Z = zeldovich(F, T, 4π/3 * ρ * a^3 / m)
        
    return 4π^2 * r^2 * a^2 * Z * (Φ / ν) * exp(exp_argument)
end