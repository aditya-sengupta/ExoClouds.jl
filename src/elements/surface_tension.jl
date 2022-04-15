function surface_tension(e::Element, T::Temperature; kwargs...)
    # if e is water, this is for water-air
    # and is equivalent to the FORTRAN "surfctwa"

    return (surface_intercept(e) - surface_slope(e) * ustrip(uconvert(K, T))) * erg/cm^2
end

# all these intercepts have been converted to a Kelvin reference
# so they may be a bit different than CARMA/setupgkern has them
surface_intercept(::water) = 118.43825 
surface_slope(::water) = -0.155

surface_intercept(::ice) = 141.0
surface_slope(::ice) = -0.15

surface_intercept(::KCl) = 179.5205
surface_slope(::KCl) = -0.07

surface_intercept(::ZnS) = 860.0
surface_slope(::ZnS) = 0.0

surface_intercept(::Na₂S) = 1033.0
surface_slope(::Na₂S) = 0.0

surface_intercept(::MnS) = 2326.0
surface_slope(::MnS) = 0.0

surface_intercept(::Cr) = 2068.63
surface_slope(::Cr) = -0.2

surface_intercept(::Fe) = 2565.2285
surface_slope(::Fe) = -0.39

surface_intercept(::Mg₂SiO₄) = 436.0
surface_slope(::Mg₂SiO₄) = 0.0

surface_intercept(::TiO₂) = 535.124
surface_slope(::TiO₂) = -0.04396

surface_intercept(::Al₂O₃) = 690.0
surface_slope(::Al₂O₃) = 0.0

# I can't dispatch this one away because it's not an air_* thing
function surface_tension_water_ice(::H₂O, T::Temperature; kwargs...)
    Tc = ustrip(uconvert(°C, T))
    return (28.5 + 0.25 * Tc) * erg/cm^2
end

function akelvin(e::Element, T::Temperature; wtpct::Float64=0.0) 
    2 * molar_weight(e) * surface_tension(e, T; wtpct=wtpct) / (R * T * density(e; wtpct=wtpct, T=T))
end


