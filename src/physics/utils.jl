
function l(η::DynamicViscosity, ρₐ::Density, μₐ::Mass, T::Temperature)
    return (2 * η) / ρₐ * sqrt(π * μₐ / (8 * R * T))
end

Kn(η::DynamicViscosity, ρₐ::Density, μₐ::Mass, T::Temperature, r::Length) = l(η, ρₐ, μₐ, T) / r

β(Kn) = 1 + 1.246*Kn + 0.42*Kn*exp(-0.87/Kn)
