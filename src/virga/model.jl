using Interpolations: LinearInterpolation
using DifferentialEquations

abstract type VfallSolution end
abstract type RGSolution end

function find_rg(f, fsed, rg, rw, α)
    # here, f = f(r, rg), and we integrate over the first argument
    fsed - moment(r -> f(r, rg), 3+α) / (rw^α) / moment(r -> f(r, rg), 3)
end

function vfall_find_root(
    ::OGVfall
    r::Length,
    atm::Atmosphere,
    T::Temperature,
    e::Element,
    w_convect::Velocity
)
    fall_velocity(r, T, atm.surface_gravity, density(e), atm.rho(atm.zref), atm.mmw, viscosity(atm, T)) - w_convect
end

@with_kw struct OGVfall <: VfallSolution
    rlo::Length{Float64}=1e-10cm,
    rhi::Length{Float64}=10.0cm
end

struct ForceBalance <: VfallSolution

function get_rw(v::OGVfall, atm::Atmosphere, T::Temperature, e::Element, w_convect::Velocity)
    find_zero(
        r -> vfall_find_root(v, r, atm, T, e, w_convect),
        (v.rlo, v.rhi),
        Roots.Brent()
    )
end

function get_rw(::ForceBalance)
    throw("not implemented")
end

struct AnalyticalRG <: RGSolution end

function calc_rg(::AnalyticalRG, atm::Atmosphere, c::VirgaCache, e::Element, p::Pressure, T::Temperature, α, qc, dz, rw, f=nothing)
    lnsig2 = (1/2) * log(c.sig)^2
    #     EQN. 13 A&M 
    #   geometric mean radius of lognormal size distribution
    rg = (c.fsed^(1.0/α) * rw * exp(-(alpha + 6) * lnsig2))

    #   droplet effective radius (cm)
    reff = rg * exp(5 * lnsig2)

    #      EQN. 14 A&M
    #   column droplet number concentration (cm^-2)
    rho_atm = p * gas_constant(atm, atm.zref) / T
    ndz = (3 * rho_atm * qc * dz /
                (4π * density(e) * rg^3) * exp(-9 * lnsig2))
    return rg, reff, ndz
end

@with_kw struct NumericalRG <: RGSolution 
    rlo::Length{Float64}=1e-10cm
    rhi::Length{Float64}=1e2cm
end

function calc_rg(n::NumericalRG, atm::Atmosphere, c::VirgaCache, e::Element, p::Pressure, T::Temperature, α, qc, dz, rw, f)
    rg = find_zero(
        rg -> find_rg(f, c.fsed, rg, rw, α),
        (n.rlo, n.rhi),
        Roots.Brent()
    )

    f_fit = (r -> f(r, rg))
    reff = moment(f_fit, 3) / moment(f_fit, 2)

    #   column droplet number concentration (cm^-2)
    rho_atm = p * gas_constant(atm, atm.zref) / T
    ndz = (3 * c.fsed * rw^α * qc * rho_atm * dz) / (4π * density(e) * moment(f_fit, 3 + α))
    return rg, reff, ndz
end


"""
dq / dz = -fsed * qc(q, z) / mixlength(q, z)
Ackerman and Marley (2001) equation 4
This treats each condensate independently, so we express that structurally
by passing in a element to this function, and broadcasting this over each
element we care about.

The rest of calc_qc.

e : the condensate being considered

Direct parameters here should only be things that could vary run to run: everything else gets passed off to VirgaCache to keep the scope clean.
"""
function optical_for_layer(
    atm::Atmosphere,
    cache::VirgaCache,
    e::Element,
    vfallsolver::VfallSolution=OGVfall(),
    rgsolver::RGSolution=AnalyticRG();
    rmin=1e-8cm, nrad=60
)
    zp = atm.zp
    nz = length(each(atm.P))
    function AM4(z::Float64, q::Float64)
        # try making this mutable for speed? 
        # i doubt it'll help in the scalar case, but a benchmark test case may be useful
        p = atm.P(z)
        T = c.temperature(z)
        qc_val = max(0.0, q - qvs(atm, T, p))
        -cache.fsed * qc_val / mixing_length(atm, T, p)
    end
    prob = ODEProblem(AM4, mixing_ratio(e, atm), (zp[1], zp[end]))
    sol = solve(prob, RK23())
    qt = LinearInterpolation(zp, sol.x, extrapolation_bc=Flat())
    mixl_out = zeros(n)
    qt_out = qt.(z)
    qc_out = maximum.(0.0, qt_out .- qvs.(atm, e, each(c.temperature), each(atm.P)))
    mixl_out = mixing_length.(atm, c, zp)
    dz = cat(1e-8cm, diff(zp), dims=1)
    rw = zeros(nz)
    rg = zeros(nz)
    reff = zeros(nz)
    ndz = zeros(nz)
    qc_path = 0.0

    for (i, z) in enumerate(zp)
        if qc_out[i] != 0.0 # layer is cloud free
            P = atm.P(z)
            T = c.temperature(z)
            K = c.Kzz(z)
            w_convect = K / mixl_out[i]
            rw[i] = get_rw.(vfallsolver, atm, T, e, w_convect)
            r_, _, _ = get_r_grid(r_min = rmin, n_radii = nrad)
            vfall_temp = get_rw_temp.(vfallsolver, r_)
            pars = linregress(r_, log.(vfall_temp)) |> coef
            α = pars[1]
            rg[i], reff[i], ndz[i] = calc_rg(rgsolver, atm, c, e, p, T, α, qc, dz, rw) # TODO handle arbitrary f
        end
        if i > 1
            qc_path = (qc_path + qc_out[i-1] *  (p_out[i-1] - p_out[i]) / atm.surface_gravity)
        end
    end

    return ([reverse(l) for l in [qc_out, qt_out, rg, reff, ndz, dz, mixl_out]])..., qc_path

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