# this has to be included after vaporpressure.jl, because it relies on the coeffs there

function latent_heat_melt(e::Element, T::Temperature) 
    latent_heat_evap(e, T)
end

latent_heat_evap(e::Element, T::Temperature) = vaporslope(e) * log(10) * R / molar_weight(e) 
# I can't find this formula in the paper referenced in the CARMA comments (Charnay et al. 2015, ApJL 813, L1), and its model also doesn't seem to use most of the condensates to which we're applying the formula. So I'll take this on faith for now but I'm skeptical (especially with the log 10).

latent_heat_evap(::H₂O, T::Temperature) = 2.501e10 * cm^2/s^2
function latent_heat_melt(::H₂O, T::Temperature) 
    Tc = ustrip(uconvert(°C, T))
    return (79.7 + 0.485 * Tc - 2.5e-3 * Tc^2) * 4.186e7 * cm^2/s^2
    # or the constant 3.337e9 cm^2/s^2
end

function latent_heat_evap(::H₂SO₄, T::Temperature; wtpctf=0.1)
    @warn "This function needs a reference to the function wtpct_tabaz, which is not in this submodule, so for now the factor wtpctf has just been arbitrarily set to 0.1. This *should* get the latent heat in the ballpark of that of water, which is also a valid approximation."
    return (1364.93 * wtpctf^3 - 1226.46 * wtpctf^2 + 382.23 * wtpctf + 540.52) * 4.184e7 * cm^2/s^2
end