using LinearAlgebra: ⋅

function l(η::DynamicViscosity, ρₐ::Density, μₐ::Mass, T::Temperature)
    return (2 * η) / ρₐ * sqrt(π * μₐ / (8 * R * T))
end

Kn(η::DynamicViscosity, ρₐ::Density, μₐ::Mass, T::Temperature, r::Length) = l(η, ρₐ, μₐ, T) / r

β(Kn) = 1 + 1.246*Kn + 0.42*Kn*exp(-0.87/Kn)

# enhancement: dispatch on sphere, hexagon, etc and make this generic
function reynolds_number(r::Length, ρₐ::Density, v::Velocity, η::DynamicViscosity)
    2 * r * ρₐ * v / η
end

function stokes_fall_velocity(r::Length, T::Temperature, g::Acceleration, ρₚ::Density, ρₐ::Density, μ::Mass, η::DynamicViscosity)
    return (2/9) * (ρₚ - ρₐ) * g * r^2 * β(Kn(η, ρₐ, μ, T, r)) / η
end

function stokes_fall_velocity(r::Length, g::Acceleration, Δρ::Density, β::Real, η::DynamicViscosity)
    return (2/9) * Δρ * g * r^2 * β / η
end

function kinetic_sedimentation_velocity(r::Length, T::Temperature, g::Acceleration, ρₚ::Density, ρₐ::Density, μₐ; A::Real=0.5)
    return A * ρₚ * g * r / (ρₐ) * sqrt(π * μₐ / (2 * R * T))
end

function fall_velocity(r::Length, T::Temperature, g::Acceleration, ρₚ::Density, ρₐ::Density, μ::Mass, η::DynamicViscosity; cdrag::Real=0.45)
    βv = β(Kn(η, ρₐ, μ, T, r))
    Δρ = ρₚ - ρₐ
    v_stokes = stokes_fall_velocity(r, g, Δρ, βv, η)
    Re = reynolds_number(r, ρₐ, g, η)
    if Re < 1
        return v_stokes
    elseif Re <= 1e3
        cd_nre2 = 32.0 * r^3 * Δρ * ρₐ * g / (3.0 * η^2) 
        b = [-0.318657e1, 0.992696, -.153193e-2, -.987059e-3, -.578878e-3, 0.855176e-4, -0.327815e-5]
        reynolds_turbulent = b ⋅ (log(cd_nre2)) .^ (0:6)
        return η * reynolds_turbulent / (2 * r * ρₐ)
    else
        return βv * sqrt(8 * Δρ * r * g / (3 * cdrag * ρₐ))
    end
end