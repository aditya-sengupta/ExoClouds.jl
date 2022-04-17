using Unitful: R
using Unitful: Length, Acceleration, Mass, Temperature, Pressure, DynamicViscosity, Density
using Unitful: ustrip
using LinearAlgebra: ⋅

using ..Elements

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

function vfall(
    r::Length,
    grav::Acceleration,
    mw_atmos::Mass,
    mfp::Length,
    visc::DynamicViscosity,
    T::Temperature,
    p::Pressure,
    rho::Density
    )

    cdrag = 0.45 

    #In order to solve the drag problem we fit y=log(reynolds)
    #as a function of x=log(cdrag * reynolds^2)
    #if you assume that at reynolds= 1, cdrag=24 and 
    #reynolds=1000, cdrag=0.45 you get the following fit: 
    # y = 0.8 * x - 0.1 * x^2
    #Full explanation: see A & M Appendix B between eq. B2 and B3
    #Simply though, this allows us to get terminal fall velocity from 
    #reynolds number
    b1 = 0.8 
    b2 = -0.01 

    R_GAS = 8.3143e7 

    #calculate constants need to get Knudsen and Reynolds numbers
    knudsen = mfp / r
    rho_atmos = p / ( (R_GAS/mw_atmos) * T )
    drho = rho - rho_atmos

    #Cunningham correction (slip factor for gas kinetic effects)
    #Cunningham, E., "On the velocity of steady fall of spherical particles through fluid medium," Proc. Roy. Soc. A 83(1910)357
    #Cunningham derived a value of 1.26 in the stone ages. In reality, this number is 
    #a function of the knudsen number. Various studies have derived 
    #different value for this number (see this citation
    #https://www.researchgate.net/publication/242470948_A_Novel_Slip_Correction_Factor_for_Spherical_Aerosol_Particles
    #Within the range of studied values, this 1.26 number changes particle sizes by a few microns
    #That is A OKAY for the level of accuracy we need. 
    beta_slip = 1. + 1.26*knudsen 

    #Stokes terminal velocity (low Reynolds number)
    #EQN B1 in A&M 
    #visc is eqn. B2 in A&M but is computed in `calc_qc`
    #also eqn 10-104 in Pruppacher & klett 1978
    vfall_r = beta_slip*(2.0/9.0)*drho*grav*r^2 / visc

    #compute reynolds number for low reynolds number case
    reynolds = 2.0*r*rho_atmos*vfall_r / visc

    #if reynolds number is between 1-1000 we are in turbulent flow 
    #limit
    if (1 < reynolds < 1e3)
        cd_nre2 = 32.0 * r^3.0 * drho * rho_atmos * grav / (3.0 * visc ^ 2 ) 
        #coefficients from EQN 10-111 in Pruppachar & Klett 1978
        #they are an empirical fit to Figure 10-9
        xx = log(cd_nre2)
        bvals = [-0.318657e1, 0.992696, -.153193e-2, -.987059e-3, -.578878e-3, 0.855176e-4, -0.327815e-5]

        y = bvals ⋅ (xx .^ (0:6))

        reynolds = exp(y)
        vfall_r = visc * reynolds / (2.0 * r * rho_atmos)
    end

    if reynolds > 1e3# 300
        # when Reynolds is greater than 1000, we can just use 
        # an asymptotic value that is independent of Reynolds number
        # Eqn. B3 from A&M 01
        vfall_r = beta_slip * sqrt( 8.0*drho*r*grav / (3.0*cdrag*rho_atmos))
    end

    return vfall_r 
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
    pvap_test = vaporpressure(m, t_test, p_test, mh)
    fx = qv_factor * pvap_test / p_test 
    return log(fx) - log(q_below)
end

"""
Root function used to find condenstation temperature. E.g. 
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
    pv = molecular_weight(m) / mmw * vaporpressure(m, t_test, p_test, mh)
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
