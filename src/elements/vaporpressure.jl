using Unitful: cm, bar, °C, K
using UnitfulAstro: dyn

vapormhfactor(::Element) = 0.0

function vaporpressure(e::Element, T::Temperature, p::Pressure=1*bar, mh::Float64=1.0)::Pressure
    if (mh != 1) && (vapormhfactor(e) == 0.0)
        throw("Warning: no M/H Dependence in vapor pressure curve for $(typeof(e))")
    end
    vaporcurve(e, T, log10(eh))
end

# for almost every species, this is true
function vaporpressure_ice(e::Element, T::Temperature, p::Pressure=1*bar, mh::Float64=1.0)::Pressure
    vaporpressure(e, T, p, mh)
end

"""
A power-law model between temperature and vapor pressure. This is the default law, but it may be overridden for individual subtypes.
"""
function vaporcurve(e::Element, T::Temperature, logmh::Float64)
    # TODO check if vaporintercept/vaporslope have physical meanings and if so rename
    10.0 ^ (vaporintercept(e) - vaporslope(e) / T - vapormhfactor(e) * logmh) * bar
end

vaporintercept(::TiO₂) = 9.5489
vaporslope(::TiO₂) = 32456.8678K # or 34602K according to "Diana email" (exocarma/carma_condensate_mod.F90)

vaporintercept(::Cr) = 7.2688
vaporslope(::Cr) = 20353.0K

vaporintercept(::ZnS) = 12.8117
vaporslope(::ZnS) = 15873.0K
vapormhfactor(::ZnS) = 1.0

function vaporcurve(::NH₃, T::Temperature, logmh::Float64)
    exp(-(86596.0 * K^2)/T^2 - (2161.0*K)/T + 10.53) * bar
end
vapormhfactor(::NH₃) = 0.0

vaporintercept(::Na₂S) = 8.5497
vaporslope(::Na₂S) = 13889.0K
vapormhfactor(::Na₂S) = 0.5

vaporintercept(::MnS) = 11.5315
vaporslope(::MnS) = 23810.0K
vapormhfactor(::MnS) = 1.0

vaporintercept(::MgSiO₃) = 11.83
vaporslope(::MgSiO₃) = 27250.0K
vapormhfactor(::MgSiO₃) = 1.0

vaporintercept(::Mg₂SiO₄) = 14.88
vaporslope(::Mg₂SiO₄) = 32488.0K
vapormhfactor(::Mg₂SiO₄) = 1.4
function vaporpressure(e::Mg₂SiO₄, T::Temperature, p::Pressure, mh::Float64)::Pressure
    vaporcurve(e, T, log10(eh)) * (p / bar)^(-0.2)
end

vaporintercept(::KCl) = 7.6106
vaporslope(::KCl) = 11382.0K

function vaporcurve(::ice, T::Temperature, logmh::Float64=0.0; do_buck::Bool=true)
    @assert T < 0°C "ice has to be cold"
    Tc = ustrip(uconvert(°C, T))
    Tk = ustrip(uconvert(K, T))
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
end

function vaporcurve(::water, T::Temperature, logmh::Float64=0.0; do_buck::Bool=true)
    @assert T > 0°C "water has to be warm"
    Tc = ustrip(uconvert(°C, T))
    Tk = ustrip(uconvert(K, T))
    if T < 1048K
        if do_buck
            return buck.BAL * exp((buck.BBL - Tc/buck.BDL)*Tc / (Tc + buck.BCL)) * dyn / cm^2
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
        return 600e6 * dyn / cm^2
    end
end

"""function vaporpressure(::ice, T::Temperature, p::Pressure=1*bar, mh::Float64=1.0)
    Tc = ustrip(uconvert(°C, T))
    # Buck 1981
    return buck.BAI * exp((buck.BBI - Tc/buck.BDI)*Tc / (Tc + buck.BCI)) * dyn / cm^2
    # Goff 1946: return 10.0 * 10^(-9.09718 * (273.16 / Tk - 1) - 3.56654 * log10(273.16 / Tk) + 0.876793 * (1 - Tk / 273.16) + log10(6.1071)) * 100 * dyn/cm^2
end"""

vaporintercept(::Fe) = 7.09
vaporslope(::Fe) = 20833.0K

function vaporcurve(::CH₄, T::Temperature, logmh::Float64=0.0)
    tcr = methane.TCRIT
    if T < tcr * K
        C = -methane.AMR * methane.AS
        B = -methane.AMR * (eethane.ALS + methane.AS + tcr)
    else
        C = -methane.AMR * methane.AL
        B = -methane.AMR * (eethane.ALV + methane.AL * tcr)
    end
    A = methane.PCRIT * tcr^(-C) * exp(-B / tcr)
    A * (T/K)^C * exp(B * K / T) * dyn/cm^2 
end

vaporintercept(::Al₂O₃) = 17.7
vaporslope(::Al₂O₃) = 45892.6K
vapormhfactor(::Al₂O₃) = 1.66

# function vaporcurve(::H₂SO₄) end