using Unitful: g, cm, u
using Unitful: Mass

function mixing_ratio(m::Molecule, mw_atmos::Mass, mh::Real, gas_mmr=nothing)
    if isnothing(gas_mmr)
        ratio = mmr_prop(m)
    else
        ratio = gas_mmr
    end
    if mh == 1
        return ratio * molecular_weight(m) / mw_atmos
    else
        throw("Alert: no M/H dependence in $(typeof(m)) routine. Consult your local theorist to determine next steps.")
    end
end

molecular_weight(::TiO₂) = 80.0u
mmr_prop(::TiO₂) = 1.69e-7
density(::TiO₂) = 4.25g/cm^3

molecular_weight(::CH₄) = 16.0u
mmr_prop(::CH₄) = 4.9e-4
density(::CH₄) = 0.49g/cm^3

molecular_weight(::NH₃) = 17.0u
mmr_prop(::NH₃) = 1.34e-4
density(::NH₃) = 0.84g/cm^3

molecular_weight(::H₂O) = 18.0u
mmr_prop(::H₂O) = 7.54e-4
density(::H₂O) = 0.93g/cm^3

molecular_weight(::Fe) = 55.845u
mmr_prop(::Fe) = 5.78e-5
density(::Fe) = 7.875g/cm^3

molecular_weight(::KCl) = 74.5u
function mixing_ratio(m::KCl, mw_atmos::Mass, mh::Real)
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
density(::KCl) = 1.99g/cm^3

molecular_weight(::MgSiO₃) = 100.4u
mixing_ratio(::MgSiO₃, mw_atmos::Mass, mh::Real) = 2.75e-3 * mh # this should be manually override-able?
density(::MgSiO₃) = 3.192g/cm^3

molecular_weight(::Mg₂SiO₄) = 140.7u
mmr_prop(::Mg₂SiO₄) = 59.36e-6
density(::Mg₂SiO₄) = 3.214g/cm^3

molecular_weight(::MnS) = 87.00u
mmr_prop(::MnS) = 6.32e-7
density(::MnS) = 4.0g/cm^3

molecular_weight(::ZnS) = 97.46u
mmr_prop(::ZnS) = 8.40e-8
density(::ZnS) = 4.04g/cm^3

molecular_weight(::Cr) = 51.996u
function mixing_ratio(m::Cr, mw_atmos::Mass, mh::Real)
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
density(::Cr) = 7.15g/cm^3

molecular_weight(::Al₂O₃) = 101.961u
mmr_prop(::Al₂O₃) = 2.51e-6
density(::Al₂O₃) = 3.987g/cm^3

molecular_weight(::Na₂S) = 78.05u
mmr_prop(::Na₂S) = 3.97e-6
density(::Na₂S) = 1.856g/cm^3