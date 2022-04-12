using Unitful
using Unitful: R

@derived_dimension MassDiffusivity 𝐋^2*𝐓^-1 true

function Kn(D::MassDiffusivity, r::Length, M::Mass, T::Temperature)
    return (3 * D / r) * sqrt(π * M / (8 * R * T))
end

function Knₜ(D, r::Length, M::Mass, T::Temperature, κ, Cₚ, μ, ρ)
    return Kn(D, r, M, T) * κ / (r * D * ρ * (Cₚ - R / (2 * μ)))
end

function D_condensate(D, r::Length, M::Mass, T::Temperature, κ, Cₚ, μ, ρ) # D' from Gao 2018 (A18); alphas force-set to 1
    Kn_val = Kn(D, r, M, T)
    λ = (1.33 * Kn_val + 0.71) / (Kn_val + 1)
    return D / (1 + λ * Kn_val)
end

function κ_condensate(D::MassDiffusivity, r::Length, M::Mass, T::Temperature, κ, Cₚ, μ, ρ) # κₐ' from Gao 2018 (A19)
    Knₜ_val = Knₜ(D, r, M, T, κ, cₚ, μ, ρ)
    λₜ = (1.33 * Knₜ_val + 0.71) / (Knₜ_val + 1)
    return κ / (1 + λₜ * Knₜ_val)
end

function dmp_dt(molecule::Molecule, D, r, M, T, κ, Cₚ, μ, ρ, S) # Jacobson 16.13 but with the numerator of Gao 2018 (A16)
    Dp = D_condensate(D, r, M, T, κ, Cₚ, μ, ρ)
    Ak = exp(2 * M * σₛ / (ρₚ * R * T * r))
    κₐ = κ_condensate(D, r, M, T, κ, Cₚ, μ, ρ)
    num = 4π * r * Dp * pₛ * (S - Ak)
    denom = ((Dp * L * pₛ * Ak) / (κₐ * T)) * (L * m / (R * T) - 1) + R * T / m
    return num / denom
end