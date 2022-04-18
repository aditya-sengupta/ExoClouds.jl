using Unitful
using Unitful: R

@derived_dimension MassDiffusivity 𝐋^2*𝐓^-1 true

function Ft_Fv(z::Length, Dp, κₐ, cₚ, is_ice::Bool)
    
    reyn_shape = re(z)
    schmidt = rmu(z) / (mw_atmos * Dp)
    prandtl = rmu(z) * cₚ / κₐ
    function _help(n)
        x = n^(1/3) * sqrt(reyn_shape)
        if is_ice
            return (x < 1) ? 1 + 0.14*x^2 : 0.86 + 0.28*x
        else
            return (x < 1) ? 1 + 0.108*x^2 : 0.78 + 0.308*x
        end
    end
    return _help(prandtl), _help(schmidt)
end

function Kn(D::MassDiffusivity, r::Length, M::Mass, T::Temperature)
    return (3 * D / r) * sqrt(π * M / (8 * R * T))
end

function Knₜ(D::MassDiffusivity, r::Length, M::Mass, T::Temperature, κ, Cₚ, μ, ρ)
    return Kn(D, r, M, T) * κ / (r * D * ρ * (Cₚ - R / (2 * μ)))
end

function D_condensate(D::MassDiffusivity, r::Length, M::Mass, T::Temperature) # D' from Gao 2018 (A18); alphas force-set to 1
    Kn_val = Kn(D, r, M, T)
    λ = (1.33 * Kn_val + 0.71) / (Kn_val + 1)
    return D / (1 + λ * Kn_val)
end

function κ_condensate(D::MassDiffusivity, r::Length, M::Mass, T::Temperature, κ, cₚ, μ, ρ) # κₐ' from Gao 2018 (A19)
    Knₜ_val = Knₜ(D, r, M, T, κ, cₚ, μ, ρ)
    λₜ = (1.33 * Knₜ_val + 0.71) / (Knₜ_val + 1)
    return κ / (1 + λₜ * Knₜ_val)
end

function dmp_dt(element::Element, atm::Atmosphere, D::MassDiffusivity, r::Length, M::Mass, T::Temperature, κ, μ, ρ, S; is_ice=true) # Jacobson 16.13 but with the numerator of Gao 2018 (A16)
    Dp = D_condensate(D, r, M, T)
    Ak = akelvin(element, T) # exp(2 * M * σₛ / (ρₚ * R * T * r))
    κₐ = κ_condensate(D, r, M, T, κ, atm.cₚ, μ, ρ)
    Ft, Fv = Ft_Fv(z, Dp, κₐ, atm.cₚ, is_ice) 
    num = 4π * r * Dp * pₛ * (S - Ak)
    denom = ((Dp * L * pₛ) / (κₐ * T)) * (L * M / (R * T) - 1) * (1 / Ft) + R * T / (M * Fv)
    return num / denom
end