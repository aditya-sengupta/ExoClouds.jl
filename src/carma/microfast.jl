# For now, this will consist of just the function signatures for every subroutine that goes into CARMA's microfast.

using Unitful
@derived_dimension SpecificParticleRate 𝐋^-3𝐓^-1 true

"""
Calculates the total amount of condensate (ice and liquid) in a particular particle type associated with a particular gas.
This is one loop iteration of CARMA/totalcondensate.F90, to be looped/mapped in a main function.

Parameters
----------
state : ...

condensate : Element

gas : Element
    The gas undergoing growth and inducing condensation.

Returns
-------
total_ice : Mass
total_liquid : Mass
"""
function condensate_quantities(state::State, condensate::Element, gas::Element)::Tuple{Mass,Mass} end

supsat_core(gas::Element, concentration::Density, T::Temperature, pvap::Pressure) = (concentration * (R / molar_weight(gas)) * T / pvap)

# reaction saturation ratio correction for type III reactions (Helling and Woitke 2006)
supsat_core(gas::Union{Na₂S,Mg₂SiO₄,Al₂O₃}, concentration::Density, T::Temperature, pvap::Pressure) = sqrt(concentration * (R / molar_weight(gas)) * T / pvap)

"""
Calculates supersaturation for a gas over liquid water/ice.
warning: if you want to vary metallicity/pressure in the vapor pressure formula, you'll have to add in some kwargs
Assumes the cloud is an ice cloud.

Because vaporpressure = vaporpressure_ice for almost every element, this could be symbolically simplified by dispatching on water separately
    but this'll do fine for now
"""
function supersaturation(gas::Element, concentration::Density, T::Temperature; 
    relative_humidity::Float64=0.0, cloud_frac::Float64=1.0, is_ice=false, kwargs...)::Tuple{Float64,Float64}
    α = relative_humidity * (1 - cloud_frac) + cloud_frac
    pvap_liquid = vaporpressure(gas, T; kwargs...)
    supsat_liquid = min(supsat_core(gas, concentration, T, pvap_liquid) - α, 0)
    if !is_ice
        return supsat_liquid
    else
        pvap_ice = vaporpressure_ice(gas, T; kwargs...)
        supsat_ice = supsat_core(gas, concentration, T, pvap_ice) - α
        supsat_ice = min(supsat_ice, (pvap_liquid - α * pvap_ice) / pvap_ice)
        return supsat_ice
    end
end