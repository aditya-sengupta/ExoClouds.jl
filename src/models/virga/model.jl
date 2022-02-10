using Unitful: R
using Dierckx: Spline1D
using Interpolations: LinearInterpolation
using DifferentialEquations

"""
dq / dz = -fsed * qc(q, z) / mixlength(q, z)
Ackerman and Marley (2001) equation 4
This treats each condensate independently, so we express that structurally
by passing in a molecule to this function, and broadcasting this over each
molecule we care about.

# TODO fix f_rain = 1, qc = q_t and run Lunine et al model for comparison
"""
function make_AM4_problem(atm::Atmosphere, zdata::Array{Length}, m::Molecule)
    # virga's calc_qc, but only the model and not the numerical solution
    # or the steps after: I'll plan to come back to this
    function AM4(z::Float64, q::Float64)
        # try making this mutable for speed? 
        # i doubt it'll help in the scalar case, but a benchmark test case may be useful
        p = atm.pressure_at_altitude(z)
        T = atm.temp_at_altitude(z)
        qc_val = max(0, q - qvs(T, p, atm))
        -atm.fsed * qc_val / mixing_length(T, p, atm)
    end
    ODEProblem(AM4, mixing_ratio(m, atm), (zdata[1], zdata[end]))
end

"""
Refine temperature pressure profile according to maximum temperature-difference between pressure layers.
why is this the description of a function called `generate_altitude`? Shouldn't it be `refine_tp_profile` or something? I guess it does both
"""
function generate_altitude(
    temp::Vector{<:Temperature}, pres::Vector{<:Pressure}, atm::Atmosphere;
    refine_tp::Bool=true, ϵ::Float64=10.0
)
    H = T -> scale_height(T, atm)
    T_p = Spline1D(pres, temp)
    Kz_p = LinearInterpolation(pres, atm.kz) # TODO check out why atm.kz is supposed to be an array?
    if refine_tp
        println("this looks complicated and unnecessary, I'll do it later")
    end
    logp = log.(pres)
    dz = -H.(T_p.(pres[1:end-1])) .* diff(logp)
    z = cumsum(dz) # may be off by one or something
    return z, Atmosphere(
        atm.gravity,
        Kz_p.(pres),
        atm.fsed,
        atm.mh,
        atm.mean_molecular_weight,
        atm.cₚ,
        atm.molecule_diameter,
        Spline1D(z, pres),
        Spline1D(z, temp),
        T_p,
        atm.supsat
    )
end