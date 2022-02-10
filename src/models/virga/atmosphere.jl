using Unitful: k, R
using ForwardDiff
using ..Molecules

include("../../utils.jl")

struct Atmosphere
    gravity::Acceleration
    kz::Float64
    fsed::Float64
    mh::Float64 # Int64?
    mean_molecular_weight::Mass
    cₚ::AbstractFloat # erg/K/g, fix typing
    molecule_diameter::Length
    pressure_at_altitude::Function
    temp_at_altitude::Function
    temp_at_pressure::Function
    supsat::Float64
end

mixing_ratio(m::Molecule, atm::Atmosphere) = mixing_ratio(m, atm.mean_molecular_weight, atm.mh)

gas_constant(atm::Atmosphere) = R / atm.mean_molecular_weight

function mean_free_path(T::Temperature, p::Pressure, atm::Atmosphere)
    n_atmos = p / (k * T)
    1 / (sqrt(2) * (n_atmos * π * atm.molecule_diameter^2))
end

function scale_height(T::Temperature, atm::Atmosphere)
    gas_constant(atm) * T / atm.grav
end

atmosphere_density(T::Temperature, p::Pressure, atm::Atmosphere) = p / (gas_constant(atm) * T)
# dHdP to be handled with autodiff
lapse_ratio(T::Temperature, p::Pressure, atm::Atmosphere) = p * derivative(atm.temp_at_pressure, p) / (T / atm.cₚ)
mixing_length(T::Temperature, p::Pressure, atm::Atmosphere) = max(0.1, lapse_ratio(T, p, atm)) * scale_height(T, atm)
qvs(T::Temperature, p::Pressure, m::Molecule, atm::Atmosphere) = (atm.supsat + 1) * vaporpressure(m, T, p) / (gas_constant(atm) * T / (atmosphere_density(T, p, atm)))