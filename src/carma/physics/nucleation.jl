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
function hom_nucleation(particle::Particle, e::Element, gas_conc::Density, T::Temperature)
    S = supersaturation(e, gas_conc, T; is_ice=particle.is_ice)
    ρ, m = density(e), molecular_weight(e)
    a = critical_radius(m, σₛ, ρ, S, T)
    σₛ = surface_tension(e, T)
    F = (4π/3) * σₛ * a^2
    Φ = diffusion_rate(e, gas_conc * R * T, T)
    cmass = 4π/3 * ρ * a^3
    Z = zeldovich(F, T, cmass / m)
    n = gas_conc / molecular_weight(e)
    return 4π * a^2 * Φ * Z * n * exp(-F / (k * T))
end

"""
The heterogeneous nucleation rate is independent of the element being nucleated onto, other than via the radius r = r(ibin,igroup). So we specify the radius as a length alongside the particle it's associated with.
"""
function het_nucleation(particle::Particle, condensing::Element, r::Length, gas_conc::Density, T::Temperature; set_zeld::Bool=false, rmiv::Float64=0.95, F_des::Energy{Float64}=2.9e-13*erg, F_sd::Energy{Float64}=2.9e-14*erg, ν::Frequency{Float64}=1e13*Hz, d_jump::Length{Float64}=1e-8*cm)
    S = supersaturation(condensing, gas_conc, T; is_ice=particle.is_ice)
    ρ, m = density(condensing), molecular_weight(condensing)
    a = critical_radius(m, σₛ, ρ, S, T)
    σₛ = surface_tension(condensing, T)
    F = (4π/3) * σₛ * a^2
    x = r / a
    ϕ = sqrt(1 - 2 * μ * x + x^2)
    f₀ = (x - μ) / ϕ
    f = (1/2) * (1 + ((1 - μ * x) / ϕ)^3 + x^3 * (2 - 3 * f₀ + f₀^3) + 3 * μ * x^2 (f₀ - 1))
    exp_argument = (2 * F_des - F_sd - F * f) / k * T
    contact_angle = acos(rmiv)
    cmass = 4π/3 * ρ * a^3

    if set_zeld
        Z = zeldovich(F, T, cmass / m)
    else
        Z = 0.1 # it's hard-coded as this in CARMA/hetnucl.F90
    end

    rnh2o = 
    
    return (Z * k * T * d_jump * rstar * sin(contact_angle) / (f * ν * m)) * 4π * r^2 * rnh2o^2 * A^2 * exp(exp_argument)
    # M is "rmw"
    # what's rstar
end