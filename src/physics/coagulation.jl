function D(η, ρₐ, μₐ, T::Temperature, r::Length)
    return k * T * β(Kn(η, ρₐ, μₐ, T, r)) / (6π * η * r)
end

function V(M::Mass, T::Temperature)
    return sqrt(8 * k * T / (π * M))
end

function λ(D, V)
    return 8 * D / (π * V)
end

"""
r - Radius
λ - Mean free path
"""
function δ(r::Length, λ::Length)
    ((2r + λ)^3 - (4r^2 + λ^2)^(3/2)) / (6*r*λ) - 2r
end

function coagulation_kernel(r₁::Length, r₂::Length, η, ρₐ, μₐ, m₁::Mass, m₂::Mass, T::Temperature)
    D₁, D₂ = D(η, ρₐ, μₐ, T, r₁),  D(η, ρₐ, μₐ, T, r₂)
    V₁, V₂ = V(m₁, T), V(m₂, T)
    λ₁, λ₂ = λ(D₁, V₁), λ(D₂, V₂)
    δ₁, δ₂ = δ(r₁, λ₁), δ(r₂, λ₂)
    num = 4π * (D₁ + D₂) * (r₁ + r₂)
    den1 = (r₁ + r₂) / (r₁ + r₂ + sqrt(δ₁^2 + δ₂^2))
    den2 = 4 * (D₁ + D₂) / ((r₁ + r₂) * sqrt(V₁^2 + V₂^2))
    return num / (den1 + den2)
end
