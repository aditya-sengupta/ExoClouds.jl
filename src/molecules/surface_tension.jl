function surface_tension(m::Molecule, T::Temperature; kwargs...)
    # if m is water, this is for water-air
    # and is equivalent to the FORTRAN "surfctwa"

    return (surface_intercept(m) - surface_slope(m) * ustrip(uconvert(K, T))) * erg/cm^2
end

# all these intercepts have been converted to a Kelvin reference
# so they may be a bit different than CARMA/setupgkern has them
surface_intercept(::H₂O) = 118.43825 
surface_slope(::H₂O) = -0.155

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

# if this comes up more, you could dispatch on ice as its own molecule, but for now it's fine
function surface_tension_water_ice(::H₂O, T::Temperature)
    Tc = ustrip(uconvert(°C, T))
    return 28.5 + 0.25 * Tc
end

function surface_tension_ice_air(::H₂O, T::Temperature)
    Tk = ustrip(uconvert(K, T))
    return 141 - 0.15 * Tk
end

function akelvin_ice(m::Molecule, T::Temperature)
    akelvin(m, T)
end

akelvin_ice(m::H₂O, T::Temperature) = 2 * molar_weight(m) * surface_tension_ice_air(m, T) / (R * T * density_ice(m))

akelvin(m::Molecule, T::Temperature; wtpct::Real=0) = 2 * molar_weight(m) * surface_tension(m, T; wtpct=wtpct) / (R * T * density(m; wtpct=wtpct, T=T))


