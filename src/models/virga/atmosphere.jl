using Unitful: cm, K # units
using Unitful: k, R, G # constants
using ForwardDiff
using Parameters

using ..Molecules

include("../../utils.jl")

gravity(mass, radius) = G * mass / (radius ^ 2)

@with_kw struct Atmosphere
    gravity::Acceleration
    kz::Real
    pressure_at_altitude::Function 
    temp_at_altitude::Function
    temp_at_pressure::Function
    # are these scaling relations always going to be splines/interpolations?
    # if so there's probably some opportunity for optimization
    # or just keeping the data around (better for sensitivity analysis)
    fsed::Real = 0.5
    mh::Real = 1.0
    mean_molecular_weight::Mass = 2.2u"u" # TODO unit check i doubt it's this lmao
    cₚ::Real = 7//2 # erg/K/g, fix typing?
    molecule_diameter::Length = 2.827e-8 * cm
    supsat::Real = 0.0
    ϵ::Real = 0.01
    ϵₖ::Temperature = 59.7 * K
end

mixing_ratio(m::Molecule, atm::Atmosphere) = mixing_ratio(m, atm.mean_molecular_weight, atm.mh)
gas_constant(atm::Atmosphere) = R / atm.mean_molecular_weight

function mean_free_path(T::Temperature, p::Pressure, atm::Atmosphere)
    n_atmos = p / (k * T)
    1 / (sqrt(2) * (n_atmos * π * atm.molecule_diameter^2))
end

scale_height(T::Temperature, atm::Atmosphere) = gas_constant(atm) * T / atm.grav
atmosphere_density(T::Temperature, p::Pressure, atm::Atmosphere) = p / (gas_constant(atm) * T)
# dHdP to be handled with autodiff
lapse_ratio(T::Temperature, p::Pressure, atm::Atmosphere) = p * derivative(atm.temp_at_pressure, p) / (T / atm.cₚ)
mixing_length(T::Temperature, p::Pressure, atm::Atmosphere) = max(0.1, lapse_ratio(T, p, atm)) * scale_height(T, atm)
qvs(T::Temperature, p::Pressure, m::Molecule, atm::Atmosphere) = (atm.supsat + 1) * vaporpressure(m, T, p) / (gas_constant(atm) * T / (atmosphere_density(T, p, atm)))