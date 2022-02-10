using Unitful: R
using Dierckx: Spline1D
using DifferentialEquations

function make_AM4_problem(atm::Atmosphere)
    """
    dq / dz = -fsed * qc(q, z) / mixlength(q, z)
    """
    function AM4!(dqdz::Float64, z::Float64, q::Float64)
        p = atm.pressure_at_altitude(z)
        T = atm.temp_at_altitude(z)
        qc_val = max(0, q - qvs(T, p, atm))
        dqdz = -atm.fsed * qc_val / mixing_length(T, p, atm)
    end
    ODEProblem(AM4, )
end

"""
Refine temperature pressure profile according to maximum temperature-difference between pressure layers.
"""
function altitude(
    temp::Vector{<:Temperature}, pres::Vector{<:Pressure}, atm::Atmosphere;
    refine_tp::Bool=true, ϵ::Float64=10.0
)
    r_atmos = R / atm.mean_molecular_weight
    H = T -> scale_height(T, r_atmos, atm.grav)

end