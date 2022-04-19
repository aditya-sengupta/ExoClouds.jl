using Unitful
using Unitful: bar, K, erg, g
using Unitful: R, k, Na, σ, c, h
using Unitful: Temperature, Mass, Length, Density, DynamicViscosity, Acceleration, Pressure, Energy, Frequency, Velocity, KinematicViscosity
using Unitful: 𝐓, 𝚯, 𝐋, 𝐌
@derived_dimension TemperatureFlux 𝚯*𝐓^-1 true
@derived_dimension SpecificParticleRate 𝐋^-3*𝐓^-1 true
@derived_dimension MassDiffusivity 𝐋^2*𝐓^-1 true
