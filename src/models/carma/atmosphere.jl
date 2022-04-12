using Unitful: G, R

@derived_dimension PressureGradient 𝐌*𝐋^-2*𝐓^-2 true

gravity(mass, radius) = G * mass / (radius ^ 2)

abstract type Horizontal end
abstract type Vertical end

struct HorizontalCartesian <: Horizontal end
xmet(::HorizontalCartesian, z::Length) = 1.0
ymet(::HorizontalCartesian, z::Length) = 1.0

struct LatitudeLongitude <: Horizontal
    planet_radius::Length
end

xmet(::LatitudeLongitude, z::Length) = planet_radius * π / 180 * cosd(yc(z))
ymet(::LatitudeLongitude, z::Length) = planet_radius * π / 180 

struct VerticalCartesian <: Vertical end
zmet(::VerticalCartesian, z::Length, dPdz::PressureGradient) = 1.0

struct Sigma <: Vertical end # Hybrid is the same as Sigma as far as I can tell
zmet(::Sigma, z::Length, dPdz::PressureGradient, g::Acceleration, ρ::Density) = abs(dPdz) / (g * ρ)

"""
The structural model aspects of a CARMA run.
"""
@with_kw struct Atmosphere # unify this with Virga's at some point
    planet_radius::Length
    surface_gravity::Acceleration
    xy::Horizontal
    z::Vertical
    mean_molecular_weight::Mass = 2.2u
    # metallicity::Real
end

gravity(atm::Atmosphere, z::Length)::Acceleration = atm.surface_gravity * atm.planet_radius^2 / (atm.planet_radius + z)^2

rhoa(atm::Atmosphere, p::Pressure, T::Temperature)::Density = p / (R / (atm.mean_molecular_weight) * T)

xmet(atm::Atmosphere, z::Length) = xmet(atm.xy, z)
ymet(atm::Atmosphere, z::Length) = ymet(atm.xy, z)
zmet(atm::Atmosphere, z::Length, p::Pressure, T::Temperature, dPdz::PressureGradient) = zmet(atm.z, z, dPdz, gravity(atm, z), rhoa(atm, p, T))
