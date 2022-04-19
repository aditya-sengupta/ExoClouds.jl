using LinearAlgebra: ⋅

using ..Elements
using Distributions: UnivariateDistribution

function advdiff_const(
    qₜ::AbstractFloat,
    ad_qbelow::AbstractFloat,
    ad_qvs::AbstractFloat,
    ad_mixl::Length,
    ad_dz::Length,
    ad_rainf::AbstractFloat,
    zb::Length,
    b::AbstractFloat,
    ϵ::AbstractFloat
)
    ad_qc = maximum(0.0, qₜ - ad_qvs)
    advdif = ad_qbelow * exp(-ad_rainf * ad_qc * ad_dz / (qₜ * ad_mixl))
    advdif - qₜ
end

function advdiff_exp(
    qₜ::AbstractFloat,
    ad_qbelow::AbstractFloat,
    ad_qvs::AbstractFloat,
    ad_mixl::Length,
    ad_dz::Length,
    ad_rainf::AbstractFloat,
    zb::Length,
    b::AbstractFloat,
    ϵ::AbstractFloat
) 
    fsed = ad_rainf
    mixl = ad_mixl
    z = ad_dz
    qc = (ad_qbelow - ad_qvs) * exp(-b * fsed / mixl * exp(zb / b) * (exp(z/b) - 1) + ϵ*z/b)
    advdif = qc + ad_qvs
    advdif - qₜ
end

vfall_find_root(r, grav, mw_atmos, mfp, visc, T, p, rho, w_convect) = vfall(r, grav, mw_atmos, mfp, visc, T, p, rho) - w_convect

""""
Define force balance for spherical particles falling in atmosphere, namely equate
gravitational and viscous drag forces.
Viscous drag assumed to be quadratic (Benchaita et al. 1983) and is a function of 
the Reynolds number dependent drag coefficient.
Drag coefficient taken from Khan-Richardson model (Richardson et al. 2002) and
is valid for 1e-2 < Re < 1e5.
Parameters
----------
vf - Quantity{<:Number, 𝐋𝐓⁻¹}
    particle sedimentation velocity (cm/s)
r - Length
    particle radius (cm)
grav - Acceleration
    acceleration of gravity (cm/s^2)
mw_atmos - Mass
    atmospheric molecular weight (g/mol)
mfp - Length
    atmospheric molecular mean free path (cm)
visc - float 
    atmospheric dynamic viscosity (dyne s/cm^2) see Eqn. B2 in A&M
t : Temperature 
    atmospheric temperature (K)
p : Quantity{<:Number,𝐌𝐋⁻¹𝐓⁻²}
    atmospheric pressure (dyne/cm^2)
rho : Quantity{<:Number,𝐌𝐋⁻³} 
    density of particle (g/cm^3)
"""
function force_balance(vf, r, grav, mw_atmos, mfp, visc, t, p, rho, gas_kinetics=true)
    rho_atmos = p / ( (R/mw_atmos) * t )
    # coefficients for drag coefficient taken from Khan-Richardson model (Richardson et al. 2002)
    # valid for 1e-2 < Re < 1e5
    a1 = 1.849; b1 = -0.31
    a2 = 0.293; b2 = 0.06
    b3 = 3.45    

    # include gas kinetic effects through slip factor 
    knudsen = mfp / r
    beta_slip = 1. + 1.26*knudsen 
    if gas_kinetics
        vf = vf/beta_slip
    end
                 
    # Reynolds number
    Re = rho_atmos * vf * 2 * r /  visc

    # Khan-Richardson approximation for drag coefficient is valid for 1e-2 < Re < 1e5
    RHS = (a1 * Re^b1 + a2 * Re^b2)^b3 * rho_atmos * π * r^2 * vf^2 
       
    # gravitational force
    LHS = 4 * π * r^3 * (rho - rho_atmos) * grav / 3 

    return LHS - RHS
end

# solve_force_balance goes into simulator

"""
Calculates the pressure of saturation mixing ratio for a gas that is extrapolated below the model domain. 
This is specifically used if the gas has saturated below the model grid
"""
function qvs_below_model(
    p_test, 
    qv_dtdlnp, 
    qv_p, 
    qv_t,
    qv_factor,
    e::Element,
    mh,
    q_below=nothing)
    #  Extrapolate temperature lapse rate to test pressure

    t_test = qv_t + log(qv_p / p_test)* qv_dtdlnp
    
    #  Compute saturation mixing ratio
    pvap_test = vaporpressure(e, t_test, p_test, mh)
    fx = qv_factor * pvap_test / p_test 
    return log(fx) - log(q_below)
end

"""
Root function used to find condensation temperature. E.g. 
the temperature when  
log p_vap = log partial pressure of gas 
Parameters
----------
t_test : float 
    Temp (K)
p_test : float 
    Pressure bars 
mh : float 
    NON log mh .. aka MH=1 for solar 
mmw : float 
    mean molecular weight (2.2 for solar)
gas_name : str 
    gas name, case sensitive 
"""
function find_cond_t(t_test::Temperature, p_test::Pressure, mh::Float64, mmw::Mass, e::Element)   
    #get vapor pressure and correct for masses of atmo and gas 
    pv = molecular_weight(e) / mmw * vaporpressure(e, t_test, p_test, mh)
    #get partial pressure
    partial_p = mixing_ratio(m, mmw, mh) * p_test * mh 
    if pv / bar == 0.0
        return -30 - log10(partial_p / bar)
    else
        return log10(pv / bar) - log10(partial_p / bar)
    end
end

function find_rg(fsed, rw, α, d::UnivariateDistribution)
    # params rg, s, loc get folded into the distribution 'd'
    fsed - moment(d, 3+α) / (rw^α) / moment(d, 3)
end
