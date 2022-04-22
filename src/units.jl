using Unitful
using Unitful: bar, K, erg, g, cm, s
using Unitful: R, k, Na, σ, c, h
using Unitful: Temperature, Mass, Length, Density, DynamicViscosity, Acceleration, Pressure, Energy, Frequency, Velocity, KinematicViscosity
using Unitful: 𝐓, 𝚯, 𝐋, 𝐌
@derived_dimension TemperatureFlux 𝚯*𝐓^-1 true
@derived_dimension SpecificParticleRate 𝐋^-3*𝐓^-1 true
@derived_dimension MassDiffusivity 𝐋^2*𝐓^-1 true

 FloatTemp = typeof(1.0K)
 FloatMass = typeof(1.0g)
 FloatLeng = typeof(1.0cm)
 FloatDens = typeof(1.0g/cm^3)
 FloatDyVi = typeof(1.0g/(cm*s))
 FloatAccl = typeof(1.0cm/s^2)
 FloatPres = typeof(1.0bar)
 FloatEner = typeof(1.0erg)
 FloatFreq = typeof(1.0/s)
 FloatVelo = typeof(1.0cm/s)
 FloatKiVi = typeof(1.0cm^2/s)
