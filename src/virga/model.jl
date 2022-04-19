using Unitful: R
# using Dierckx: Spline1D
using Interpolations: LinearInterpolation
using DifferentialEquations

abstract type fsedParam end

struct fsedConst <: fsedParam end
@with_kw struct fsedExpon <: fsedParam 
    ϵ::Float64 = 0.01
    b::Float64 = 1.0
end

get_fsed_mid(::fsedConst, c::Cache, fsed) = fsed
get_fsed_mid(f::fsedExpon, c::Cache, fsed) = fsed / exp(cache.z_alpha / f.b)

function dq(::fsedConst, cache::VirgaCache, fsed, q_below, q_sat, mixlength, z_bot)
    return q_sat + (q_below - q_sat) * exp.(-fsed .* cache.layer_dz ./ mixlength)
end

function dq(f::fsedExpon, cache::VirgaCache, fsed, q_below, q_sat, mixlength, z_bot)
    return q_sat + (q_below - q_sat) * exp((-b * fsed_mid(f, cache, fsed) * exp((z_bot .+ cache.layer_dz) / b - 1) + f.ϵ .* cache.layer_dz) ./ mixlength)
end

abstract type VfallSolution end

@with_kw struct OGVfall <: VfallSolution
    rlo::Length{Float64}=1e-10cm,
    rhi::Length{Float64}=10.0cm
end
struct ForceBalance <: VfallSolution

function get_rw_temp(v::OGVfall)
    function vfall_find_root(...) end
    find_zero(
        vfall_find_root,
        (v.rlo, v.rhi),
        Roots.Brent()
    )
end

function get_rw_temp(::ForceBalance)
    # solve_force_balance
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
    T_extr = LinearInterpolation((atm.zp,), Tp, extrapolation_bc=Flat())
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

Direct parameters here should only be things that could vary run to run: everything else gets passed off to VirgaCache to keep the scope clean.
"""
function optical_for_layer(
    atm::Atmosphere,
    cache::VirgaCache,
    e::Element,
    qt_below::Float64,
    fsed::Float64,
    mixlength::Length{Float64},
    z_bot::Length{Float64},
    rho_p::Density{Float64},
    f::fsedParam=fsedConst(),
    vfallsolver::VfallSolution=OGVfall()
)

    qvs_val = qvs(atm, e, T, p)
    if qt_below < qvs_val
        return 0.0 # or whatever
    end

    qt_top = dq(f, cache, fsed, qt_below, qvs_val, mixlength, z_bot)
    
    qt_layer = (1/2) * (qt_below .+ qt_top)
    qc_layer = maximum.(0, qt_layer - qvs_val)

    rw_temp = get_rw_temp(vfallsolver)
    lnsig2 = 0.5 * log(cache.sig)^2
    sig_alpha = max(1.1, sig)
    r_, _, _ = get_r_grid(...) # TODO 
    vfall_temp = zeros(length(r_))
    for j = 1:length(r_)
        vfall_temp[j] = get_vfall_temp(::vfallsolver, ...) # TODO
    end
    pars = linregress(r_, log.(vfall_temp)) |> coef
    alpha = pars[1]
    fsed_mid = get_fsed_mid(f, cache, fsed)
    rg_layer = (fsed_mid^(1/alpha)) * rw_layer * exp(-(alpha+6) * lnsig2)
    reff_layer = rg_layer * exp(5 * lnsig2)
    ndz_layer = (3*atm.rho(atm.zref) .* qc_layer .* dz_layer ./ (4π * rho_p*rg_layer^3 ) * exp(-9*lnsig2))
    return ndz_layer
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
    T_p = LinearInterpolation((pres,), temp, extrapolation_bc=Flat())
    Kz_p = LinearInterpolation((pres,), atm.kz, extrapolation_bc=Flat())
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