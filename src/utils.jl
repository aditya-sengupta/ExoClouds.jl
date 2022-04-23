using Parameters
using Unitful
using Unitful: bar, K, erg, g, cm, s, u
using Unitful: R, k, Na, σ, c, h, G
using Unitful: Temperature, Mass, Length, Density, DynamicViscosity, Acceleration, Pressure, Energy, Frequency, Velocity, KinematicViscosity
using Unitful: 𝐓, 𝚯, 𝐋, 𝐌
@derived_dimension TemperatureFlux 𝚯*𝐓^-1 true
@derived_dimension SpecificParticleRate 𝐋^-3*𝐓^-1 true
@derived_dimension MassDiffusivity 𝐋^2*𝐓^-1 true

# forces unit conversions to CGS every time a parameter is passed into a function or struct
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

using Interpolations
using Interpolations: LinearInterpolation, Extrapolation

xr = [1.0, 2.0] # dummy values

LMExtrapolation = typeof(LinearInterpolation((xr .* cm,), xr .* g, extrapolation_bc=Flat()))
LPExtrapolation = typeof(LinearInterpolation((xr .* cm,), xr .* bar, extrapolation_bc=Flat()))
LDExtrapolation = typeof(LinearInterpolation((xr .* cm,), xr .* g/cm^3, extrapolation_bc=Flat()))
LTExtrapolation = typeof(LinearInterpolation((xr .* cm,), xr .* K, extrapolation_bc=Flat()))
