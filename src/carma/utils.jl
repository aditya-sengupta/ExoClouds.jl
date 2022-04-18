using Parameters
using Unitful: Length, Acceleration, Mass, Pressure, Temperature
using Unitful: σ, k, c, h, R
using Interpolations: LinearInterpolation, Extrapolation

@derived_dimension TemperatureFlux 𝚯*𝐓^-1 true
@derived_dimension SpecificParticleRate 𝐋^-3𝐓^-1 true

function l(η, ρₐ, μₐ, T::Temperature)
    return (2 * η) / ρₐ * sqrt(π * μₐ / (8 * R * T))
end

Kn(η, ρₐ, μₐ, T::Temperature, r::Length) = l(η, ρₐ, μₐ, T) / r

β(Kn) = 1 + 1.246*Kn + 0.42*Kn*exp(-0.87/Kn)
