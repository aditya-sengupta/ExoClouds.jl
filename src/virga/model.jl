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
    f::fsedParam=fsedConst(),
    vfallsolver::VfallSolution=OGVfall()
)
    qc_path = 0.0
    qvs_val = qvs(atm, e, T, p)
    if qt_below < qvs_val
        throw("need to determine the output signature and put zeros in all of them")
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
    ndz_layer = (3*atm.rho(atm.zref) .* qc_layer .* dz_layer ./ (4π * density(e) * rg_layer^3 ) * exp(-9*lnsig2))
    return qt_top, qc_sub, qt_sub, reff_layer, ndz_layer
end

"""
Refine temperature pressure profile according to maximum temperature-difference between pressure layers, and make a Virga atmosphere object.
"""
function atmosphere_virga(
    temperature::Vector{Temperature{Float64}}, pressure::Vector{Pressure{Float64}}, Kzzs::Vector{KinematicViscosity{Float64}} mwp::Mass{Float64}, planet_radius::Length{Float64}, 
    surface_gravity::Acceleration{Float64};
    refine_tp::Bool=true, ϵ::Temperature{Float64}=10.0K,
    kwargs...
)
    n = length(pressure)
    T_p = LinearInterpolation((pressure,), temperature, extrapolation_bc=Flat())
    Kz_p = LinearInterpolation((pressure,), Kzzs, extrapolation_bc=Flat())

    if refine_tp
        temps = T_p.(pressure)
        dtemps = abs(diff(temps))
        while maximum(dtemps) > ϵ
            indx = findfirst(dtemps .> ϵ)
            mids = (pressure[indx+1] + pressure[indx]) / 2
            insert!(pressure, indx+1, mids)
        end
    end

    T = zeros(n)
    K = zeros(n)
    T[1:end-1] = T_p.(pressure[1:end-1])
    T[end] = temperature[end]
    K[1:end-1] = Kz_p.(pressure[1:end-1])
    K[end] = Kzzs[end]

    logp = log.(reverse(pres_))
    dz = -scale_height.(atm, T[1:end-1]) .* diff(logp)
    z = cumsum(dz)
    prepend!(0m, z)
    
    atm = Atmosphere(
        planet_radius,
        surface_gravity,
        z, pressure, mwp .* ones(length(z)); kwargs...
    )
    return atm, T, K
end