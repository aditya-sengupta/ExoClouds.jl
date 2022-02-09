using Unitful, UnitfulAstro
using Unitful: cm, bar, °C, K
using UnitfulAstro: dyn

include("constants.jl")
include("molecules.jl")
include("../utils.jl")

vapormhfactor(::Molecule) = 0.0

function vaporpressure(m::Molecule, T::Temperature, mh::AbstractFloat=1)
    if (mh != 1) && (vapormhfactor(m) == 0.0)
        throw("Warning: no M/H Dependence in vapor pressure curve for $(typeof(m))")
    end
    vaporcurve(m, T, log10(mh))
end

"""
A power-law model between temperature and vapor pressure. This is the default law, but it may be overridden for individual subtypes.
"""
function vaporcurve(m::Molecule, T::Temperature, logmh::AbstractFloat)
    # TODO check if vaporintercept/vaporslope have physical meanings and if so rename
    10.0 ^ (vaporintercept(m) - vaporslope(m) / T - vapormhfactor(m) * logmh) * bar
end

vaporintercept(::TiO₂) = 9.5489
vaporslope(::TiO₂) = 32456.8678u"K"

vaporintercept(::Cr) = 7.2688
vaporslope(::Cr) = 20353.0u"K"

vaporintercept(::ZnS) = 12.8117
vaporslope(::ZnS) = 15873.0u"K"
vapormhfactor(::ZnS) = 1.0

function vaporcurve(::NH₃, T::AbstractFloat, logmh::AbstractFloat)
    exp(-(86596.0 * K^2)/T^2 - (2161.0*K)/T + 10.53) * bar
end
vapormhfactor(::NH₃) = 0.0

vaporintercept(::Na₂S) = 8.5497
vaporslope(::Na₂S) = 13889.0u"K"
vapormhfactor(::Na₂S) = 0.5

vaporintercept(::MnS) = 11.5315
vaporslope(::MnS) = 23810.0u"K"
vapormhfactor(::MnS) = 1.0

vaporintercept(::MgSiO₃) = 11.83
vaporslope(::MgSiO₃) = 27250.0u"K"
vapormhfactor(::MgSiO₃) = 1.0

vaporintercept(::Mg₂SiO₄) = 14.88
vaporslope(::Mg₂SiO₄) = 32488.0u"K"
vapormhfactor(::Mg₂SiO₄) = 1.4
function vaporpressure(m::Mg₂SiO₄, T::AbstractFloat, logmh::AbstractFloat, p::AbstractFloat)
    # this one doesn't quite fit the signature of the parent, but dispatch will take care of the case where you don't want to specify a p
    vaporcurve(m, T, logmh) * p / 10 ^ (6.2)
end

vaporintercept(::KCl) = 7.6106
vaporslope(::KCl) = 11382.0u"K"

function vaporcurve(::H₂O, T::Temperature, logmh=0.0, do_buck::Bool=true)
    Tc = ustrip(uconvert(°C, T))
    Tk = ustrip(uconvert(K, T))
    if T < 0°C
        if do_buck
            # bit annoying to type, but keeps the namespace clean
            return buck.BAI * exp((buck.BBI - Tc/buck.BDI)*Tc / (Tc + buck.BCI)) * dyn / cm^2
        else
            return 10.0 * exp(
                1.0 / Tk * (
                    wexler.HH0 + (
                        wexler.HH1 + wexler.HH5 * log(Tk) + (
                            wexler.HH2 + (
                                wexler.HH3 + wexler.HH4 * Tk
                            ) * Tk
                        ) * Tk
                    ) * Tk
                )
            ) * dyn / cm^2
        end
    elseif T < 1048u"K"
        if do_buck
            return buck.BAL * exp((buck.BBL - Tc/buck.BDL)*Tc / (Tc + buck.BCL))
        else
            return 10.0 * exp(
                (1.0 * (Tk^2)) * (
                    wexler.GG0 + (
                        wexler.GG1 + (
                            wexler.GG2 + wexler.GG7 * log(Tk) + (
                                wexler.GG3 + (
                                    wexler.GG4 + (
                                        wexler.GG5 + wexler.GG6 * Tk
                                    ) * Tk
                                ) * Tk
                            ) * Tk
                        ) * Tk
                    ) * Tk
                )
            ) * dyn/cm^2
        end
    else
        return 600 * bar
    end
end

vaporintercept(::Fe) = 7.09
vaporslope(::Fe) = 20833.0u"K"

function vaporcurve(::Fe, T::Temperature, logmh=0.0)
    tcr = methane.TCRIT
    if T < tcr * K
        C = -methane.AMR * methane.AS
        B = -methane.AMR * (methane.ALS + methane.AS + tcr)
    else
        C = -methane.AMR * methane.AL
        B = -methane.AMR * (methane.ALV + methane.AL * tcr)
    end
    A = methane.PCRIT * tcr^(-C) * exp(-B / tcr)
    A * T^C * exp(B / T) # here be unit issues 
end

vaporintercept(::Al₂O₃) = 17.7
vaporslope(::Al₂O₃) = 45892.6u"K"
vapormhfactor(::Al₂O₃) = 1.66