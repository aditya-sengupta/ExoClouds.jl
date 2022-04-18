using Unitful: R
# using Dierckx: Spline1D
using Interpolations: LinearInterpolation
using DifferentialEquations

abstract type fsedParam end

struct fsedConst <: fsedParam end
struct fsedExpon <: fsedParam 
    b::Float64
end

function dq(::fsedConst, fsed, q_below, q_sat, fsed, Δz, mixlength, ϵ)
    return q_sat + (q_below - q_sat) * exp(-fsed * Δz / mixlength)
end

function dq(::fsedExpon, fsed, q_below, q_sat, fsed, Δz, mixlength, ϵ)
    fs = fsed / exp()
end

"""
dq / dz = -fsed * qc(q, z) / mixlength(q, z)
Ackerman and Marley (2001) equation 4
This treats each condensate independently, so we express that structurally
by passing in a element to this function, and broadcasting this over each
element we care about.
"""
function make_AM4_problem(atm::Atmosphere, e::Element, Tp::Vector{Temperature{Float64}})
    # virga's calc_qc, but only the model and not the numerical solution
    # or the steps after: I'll plan to come back to this
    T_extr = extrapolate(atm.zp, Tp)
    function AM4(z::Float64, q::Float64)
        # try making this mutable for speed? 
        # i doubt it'll help in the scalar case, but a benchmark test case may be useful
        p = atm.P(z)
        T = T_extr(z)
        qc_val = max(0.0, q - qvs(atm, T, p))
        -atm.fsed * qc_val / mixing_length(atm, T, p)
    end
    ODEProblem(AM4, mixing_ratio(e, atm), (atm.zp[1], atm.zp[end]))
end

"""
The rest of calc_qc.

e : the condensate being considered
qt_below : total mixing ratio below the layer
T : temperature at the midpoint of the layer
p : pressure at the midpoint of the layer
"""
function optical_for_layer(
    e::Element,
    qt_below::Float64,
    T::Temperature{Float64},
    p::Pressure{Float64}
)

    qvs_val = qvs(atm, e, T, p)
    if qt_below < qvs_val
        return 0.0 # or whatever
    end

    q_range = (qt_below / 1000, qt_below)

    
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
    T_p = extrapolate(pres, temp)
    Kz_p = extrapolate(pres, atm.kz) # TODO check out why atm.kz is supposed to be an array?
    if refine_tp
        println("refine_tp looks complicated and unnecessary, I'll do it later")
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