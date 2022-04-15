using Unitful: G, R

@derived_dimension PressureGradient 𝐌*𝐋^-2*𝐓^-2 true

gravity(mass::Mass, radius::Length) = G * mass / (radius ^ 2)

abstract type Horizontal{N} end
abstract type Vertical{N} end

struct HorizontalCartesian{N} <: Horizontal 
    ref_x::SVector{N,Float64}
    ref_y::SVector{N,Float64}
end

xmet(::HorizontalCartesian, z::Length) = 1.0
ymet(::HorizontalCartesian, z::Length) = 1.0

struct LatitudeLongitude{N} <: Horizontal
    planet_radius::Length
    ref_lat::SVector{N,Float64}
    ref_lon::SVector{N,Float64}
end

xmet(::LatitudeLongitude, z::Length) = planet_radius * π / 180 * cosd(yc(z))
ymet(::LatitudeLongitude, z::Length) = planet_radius * π / 180 

struct VerticalCartesian{N} <: Vertical 
    ref_z::SVector{N,Float64}
end
zmet(::VerticalCartesian, z::Length, dPdz::PressureGradient) = 1.0

struct Sigma{N} <: Vertical 
    ref_z::SVector{N,Float64}
end # Hybrid is the same as Sigma as far as I can tell
zmet(::Sigma, z::Length, dPdz::PressureGradient, g::Acceleration, ρ::Density) = abs(dPdz) / (g * ρ)

gravity(atm::Atmosphere, z::Length)::Acceleration = atm.surface_gravity * atm.planet_radius^2 / (atm.planet_radius + z)^2

rhoa(atm::Atmosphere, p::Pressure, T::Temperature)::Density = p / (R / (atm.mean_molecular_weight) * T)

xmet(atm::Atmosphere, z::Length) = xmet(atm.xy, z)
ymet(atm::Atmosphere, z::Length) = ymet(atm.xy, z)
zmet(atm::Atmosphere, z::Length, p::Pressure, T::Temperature, dPdz::PressureGradient) = zmet(atm.z, z, dPdz, gravity(atm, z), rhoa(atm, p, T))
