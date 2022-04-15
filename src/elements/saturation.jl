using Unitful

supsat_core(gas::Element, concentration::Density, T::Temperature, pvap::Pressure) = (concentration * (R / molar_weight(gas)) * T / pvap)

# reaction saturation ratio correction for type III reactions (Helling and Woitke 2006)
supsat_core(gas::Union{Na₂S,Mg₂SiO₄,Al₂O₃}, concentration::Density, T::Temperature, pvap::Pressure) = sqrt(concentration * (R / molar_weight(gas)) * T / pvap)

"""
Calculates supersaturation for a gas over liquid water/ice.
warning: if you want to vary metallicity/pressure in the vapor pressure formula, you'll have to add in some kwargs

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