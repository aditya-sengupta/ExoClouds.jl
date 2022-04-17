
function stokes_fall_velocity(ρₚ::Density, η, ρₐ, μₐ, g::Acceleration, r::Length)
    return (2/9) * ρₚ * g * r^2 * β(Kn(η, ρₐ, μₐ, T, r)) / η
end

function kinetic_sedimentation_velocity(ρₚ::Density, g::Acceleration, r::Length, ρₐ::Density, μₐ, T::Temperature; A=0.5)
    return A * ρₚ * g * r / (ρₐ) * sqrt(π * μₐ / (2 * R * T))
end