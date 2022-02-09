function mmr(m::Molecule, mw_atmos::Mass, mh::AbstractFloat, gas_mmr=nothing)
    if isnothing(gas_mmr)
        ratio = mmr_ratio(m)
    else
        ratio = gas_mmr
    end
    if mh == 1
        mmr_ratio(m) * mw(m) / mw_atmos
    else
        throw("Alert: no M/H dependence in $(typeof(m)) routine. Consult your local theorist to determine next steps.")
    end
end

mw(::TiO₂) = 80.0u"u"
mmr_ratio(::TiO₂) = 1.69e-7
ρ(::TiO₂) = 4.25*g/cm^3

mw(::CH₄) = 16.0u"u"
mmr_ratio(::CH₄) = 4.9e-4
ρ(::CH₄) = 0.49*g/cm^3

mw(::NH₃) = 17.0u"u"
mmr_ratio(::CH₄) = 1.34e-4
ρ(::CH₄) = 0.84*g/cm^3

mw(::H₂O) = 18.0u"u"
mmr_ratio(::H₂O) = 7.54e-4
ρ(::H₂O) = 0.93*g/cm^3

mw(::Fe) = 55.845u"u"
mmr_ratio(::Fe) = 5.78e-5
ρ(::Fe) = 7.875*g/cm^3

mw(::KCl) = 74.5u"u"
function mmr(m::KCl, mw_atmos::Mass, mh::AbstractFloat)
    if mh == 1
        ratio = 2.2627e-7
    elseif mh == 10
        ratio = 2.1829e-6
    elseif mh == 50
        ratio = 8.1164e-6
    else
        throw("KCl gas properties can only be computed for 1, 10, and 50x solar metallicity")
    end
    ratio * mw(m) / mw_atmos
end
ρ(::KCl) = 1.99*g/cm^3

mw(::MgSiO₃) = 100.4u"u"
mmr_ratio(::MgSiO₃) = 2.75e-3
ρ(::MgSiO₃) = 3.192*g/cm^3

mw(::Mg₂SiO₄) = 140.7u"u"
mmr_ratio(::Mg₂SiO₄) = 59.36e-6
ρ(::Mg₂SiO₄) = 3.214*g/cm^3

mw(::MnS) = 87.00u"u"
mmr_ratio(::MnS) = 6.32e-7
ρ(::MnS) = 4.0*g/cm^3

mw(::ZnS) = 97.46u"u"
mmr_ratio(::ZnS) = 8.40e-8
ρ(::ZnS) = 4.04*g/cm^3

mw(::Cr) = 51.996u"u"
function mmr(m::KCl, mw_atmos::Mass, mh::AbstractFloat)
    if mh == 1
        ratio = 8.87e-7
    elseif mh == 10
        ratio = 8.6803e-6
    elseif mh == 50
        ratio = 4.1308e-5
    else
        throw("Cr gas properties can only be computed for 1, 10, and 50x solar metallicity")
    end
    ratio * mw(m) / mw_atmos
end
ρ(::Cr) = 7.15*g/cm^3

mw(::Al₂O₃) = 101.961u"u"
mmr_ratio(::Al₂O₃) = 2.51e-6
ρ(::Al₂O₃) = 3.987*g/cm^3

mw(::Na₂S) = 78.05u"u"
mmr_ratio(::Na₂S) = 3.97e-6
ρ(::Na₂S) = 1.856*g/cm^3