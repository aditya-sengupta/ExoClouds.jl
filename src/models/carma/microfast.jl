# For now, this will consist of just the function signatures for every subroutine that goes into CARMA's microfast.

using Unitful
@derived_dimension SpecificParticleRate 𝐋^-3𝐓^-1 true

"""
Calculates the total amount of condensate (ice and liquid) in a particular particle type associated with a particular gas.
This is one loop iteration of CARMA/totalcondensate.F90, to be looped/mapped in a main function.

Parameters
----------
state : ...

condensate : Molecule

gas : Molecule
    The gas undergoing growth and inducing condensation.

Returns
-------
total_ice : Mass
total_liquid : Mass
"""
function condensate_quantities(state::State, condensate::Molecule, gas::Molecule)::Tuple{Mass,Mass} end

supsat_core(gas::Molecule, concentration::Density, T::Temperature, pvap::Pressure) = (concentration * (R / molar_weight(gas)) * T / pvap)

# reaction saturation ratio correction for type III reactions (Helling and Woitke 2006)
supsat_core(gas::Union{Na₂S,Mg₂SiO₄,Al₂O₃}, concentration::Density, T::Temperature, pvap::Pressure) = sqrt(concentration * (R / molar_weight(gas)) * T / pvap)

"""
Calculates supersaturation for a gas over liquid water/ice.
warning: if you want to vary metallicity/pressure in the vapor pressure formula, you'll have to add in some kwargs
Assumes the cloud is an ice cloud.
"""
function supersaturations(gas::Molecule, concentration::Density, T::Temperature; relative_humidity::Real=0.0, cloud_frac::Real=1.0)::Tuple{Real,Real}
    α = relative_humidity * (1 - cloud_frac) + cloud_frac
    pvap_liquid = vaporpressure(gas, T)
    pvap_ice = vaporpressure_ice(gas, T)
    supsat_liquid = min(supsat_core(gas, concentration, T, vaporpressure(gas, T)) - α, 0)
    supsat_ice = supsat_core(gas, concentration, T, pvap_ice) - α
    supsat_ice = min(supsat_ice, (pvap_liquid - α * pvap_ice) / pvap_ice)
    return supsat_liquid, supsat_ice
end
