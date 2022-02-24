function mixing_ratio(m::Molecule, mw_atmos::Mass, mh::AbstractFloat, gas_mmr=nothing)
    if isnothing(gas_mmr)
        ratio = mmr_prop(m)
    else
        ratio = gas_mmr
    end
    if mh == 1
        return mmr_prop(m) * molecular_weight(m) / mw_atmos
    else
        throw("Alert: no M/H dependence in $(typeof(m)) routine. Consult your local theorist to determine next steps.")
    end
end

molecular_weight(::TiO₂) = 80.0u"u"
mmr_prop(::TiO₂) = 1.69e-7
density(::TiO₂) = 4.25*g/cm^3

molecular_weight(::CH₄) = 16.0u"u"
mmr_prop(::CH₄) = 4.9e-4
density(::CH₄) = 0.49*g/cm^3

molecular_weight(::NH₃) = 17.0u"u"
mmr_prop(::NH₃) = 1.34e-4
density(::NH₃) = 0.84*g/cm^3

molecular_weight(::H₂O) = 18.0u"u"
mmr_prop(::H₂O) = 7.54e-4
density(::H₂O) = 0.93*g/cm^3

molecular_weight(::Fe) = 55.845u"u"
mmr_prop(::Fe) = 5.78e-5
density(::Fe) = 7.875*g/cm^3

molecular_weight(::KCl) = 74.5u"u"
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
    ratio * molecular_weight(m) / mw_atmos
end
density(::KCl) = 1.99*g/cm^3

molecular_weight(::MgSiO₃) = 100.4u"u"
mmr_prop(::MgSiO₃) = 2.75e-3
density(::MgSiO₃) = 3.192*g/cm^3

molecular_weight(::Mg₂SiO₄) = 140.7u"u"
mmr_prop(::Mg₂SiO₄) = 59.36e-6
density(::Mg₂SiO₄) = 3.214*g/cm^3

molecular_weight(::MnS) = 87.00u"u"
mmr_prop(::MnS) = 6.32e-7
density(::MnS) = 4.0*g/cm^3

molecular_weight(::ZnS) = 97.46u"u"
mmr_prop(::ZnS) = 8.40e-8
density(::ZnS) = 4.04*g/cm^3

molecular_weight(::Cr) = 51.996u"u"
function mmr(m::Cr, mw_atmos::Mass, mh::AbstractFloat)
    if mh == 1
        ratio = 8.87e-7
    elseif mh == 10
        ratio = 8.6803e-6
    elseif mh == 50
        ratio = 4.1308e-5
    else
        throw("Cr gas properties can only be computed for 1, 10, and 50x solar metallicity")
    end
    ratio * molecular_weight(m) / mw_atmos
end
density(::Cr) = 7.15*g/cm^3

molecular_weight(::Al₂O₃) = 101.961u"u"
mmr_prop(::Al₂O₃) = 2.51e-6
density(::Al₂O₃) = 3.987*g/cm^3

molecular_weight(::Na₂S) = 78.05u"u"
mmr_prop(::Na₂S) = 3.97e-6
density(::Na₂S) = 1.856*g/cm^3