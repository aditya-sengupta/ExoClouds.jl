using Unitful: σ, k, c, h

"""
Calculates average Planck intensity in a wavelength band.
This breaks the modeling/simulation barrier, 
so ideally I'll come back to this, find the paper, and reformulate it as a smaller diffeq.
"""
function planckBandIntensityConley2011(λ_center::Length, λ_width::Length, T::Temperature, num_partitions=1)
    total = 0

    kt = k * T
    λ_lower = λ_center - λ_width / 2
    λ_upper = λ_center + λ_width / 2

    if (num_partitions > 1)
        dellam = (λ_lower / λ_upper) * exp(1/num_partitions)
        f₁ = c / (λ_lower / dellam)
        f₂ = c / (λ_lower)
    else
        dellam = 1
        f₁ = c / λ_upper
        f₂ = c / λ_lower
    end

    for i = 1:num_partitions
        Δf = 1/2 * (f₂ - f₁)
        f_cen = 1/2 * (f₁ + f₂)
        mi = 1 / f_cen
        α = h / kt

        argexp = α * f_cen
        if (argexp < 0.5)
            e = 1 + argexp + 1/2 * argexp^2 + 1/6 * argexp^3 + 1/24 * argexp^4
            em1i = 1 / (e - 1)
            di = e * em1i
        elseif (argexp < 20)
            e = exp(argexp)
            em1i = 1 / (e - 1)
            di = e * em1i
        else
            e = 1e20 # frequency >> temperature
            em1i = 1e-20
            di = 1
        end

        # frontpiece
        coeff = 2h * f_cen^3 * em1i / c^2

        # integrals
        o = f₂ - f₁     # int 1 deps
        tt = 2/3 * Δf^3 # int eps^2 deps

        # term and 4th order correction
        t1 = 1
        t3 = 3*mi^2 - 3*α*di*mi + α^2*di^2 - 1/2 * α^2 * di

        total += coeff * (o * t1 + tt * t3)

        f₂ = f₁
        f₁ *= dellam
    end

    return total / λ_width
end

"""
Add in swelling due to humidity (currently not implemented.)
"""
function getwetr(r::Length, T::Temperature)
    return r
end

"""
CARMA's "gro": growth due to diffusion processes
"""
function growth_diffusion(m::Molecule, r::Length, T::Temperature, Dp::MassDiffusivity, fv::Real)
    br = getwetr(r, T) # low
    return 4π * br * Dp * fv * molar_weight(m) / (k * T * Na)
end

"""
CARMA's "gro1": growth due to conduction
"""
function growth_conduction(particle::Particle, m::Molecule, r::Length, T::Temperature, κₐ, ft::Real)
    rlh = latent_heat(particle, m, T)
    br = getwetr(r, T) # low 
    return molar_weight(m) * rlh^2 / (R * T^2 * ft * κₐ) / (4π * br)
end

"""
CARMA's "gro2": growth due to radiation
"""
function growth_radiation(particle::Particle, m::Molecule, T::Temperature)
    return 1 / latent_heat(particle, m, T)
end

"""
Particle loss rate due to heating.

rad_incident is FORTRAN's radint, to be set by the user.

This should be called at a particular (r, T) which the state knows.
"""
function pheat(atm::Atmosphere, particle::Particle, m::Molecule, r::Length, T::Temperature, p::Pressure, conc::Density, relative_humidity::Real, cloud_frac::Real, rad_incident=0, λ::Vector=[], dλ::Vector=[];
    NWAVE=0, POWMAX=85, wtpct=0.1, do_pheat=true, do_pheatatm=true, do_mie=true, do_wave_emit=true, max_iter=10) 
    if particle.is_ice
        expon = max(-POWMAX, akelvin_ice(m, T) / getwetr(r, T))
        akas = exp(expon)
    elseif (typeof(m) == H₂SO₄)
        @warn "make sure wtpct is set correctly"
        argsol = wtpct / 100 / molar_weight(m) / ((1 - wtpct / 100) / molar_weight(water) + wtpct/100/molar_weight(m))
        expon = max(-POWMAX, akelvin(m, T) / getwetr(r, T))
        akas = exp(expon) * argsol
    else
        expon = max(-POWMAX, akelvin(m, T) / getwetr(r, T))
        akas = exp(expon)
    end

    # Dp, κₐ, fv, ft?
    g0 = growth_diffusion(m, r, T, Dp, fv)
    g1 = growth_conduction(particle, m, r, T, κₐ, ft)
    g2 = growth_radiation(particle, m, T)

    pvap = vaporpressure(particle, m, T, p, atm.mh)
    ss = supersaturations(m, conc, T; relative_humidity=relative_humidity, cloud_frac=cloud_frac)[particle.is_ice ? 2 : 1]
    
    if (!do_pheat || !do_mie)
        dm_dt = pvap * (ss + 1 - akas) * g0 / (1 + g0 * g1 * pvap)
    else
        rlh = latent_heat(particle, m, T)
        dtp = 0.0 # I've taken out some limiting
        qrad = 0.0
        for iter = 1:max_iter
            for iwvl = 1:NWAVE
                plkint = do_wave_emit ? planckBandIntensityConley2011(λ[iwvl], dλ[iwvl], T) : 0.0
                qrad += 4π * (1 - ssa) * qext * π * r^2 * (rad_incident - plkint * dλ[iwvl])
            end
            dm_dt = pvap * (ss + 1 - akas * (1 + qrad * g1 * g2)) * g0 / (1 + g0 * g1 * pvap)
            # calculate new particle temperature
            if dm_dt < -drmass_dt # ibin + 1
                dtp = ((rlh * dm_dt) + qrad) / (4π * getwetr(r, T) * κₐ * Ft)
            end
        end
        if do_pheatatm
            if dm_dt > -drmass_dt
                phprod += 4π * getwetr(r, T) * κₐ * Ft * dtp * conc / (atm.cₚ * ρ_atm)
            end
        end
    end
    return dtp
end

"""
FORTRAN loops this over the particle groups, which we'll do in the main loop.
This is a function of just "particle", which translates in FORTRAN to "igroup".
"""
function growevapl(state::State, particle::Molecule)
    # can't finish this yet because ratt1,2,3 aren't defined as far as I can see
    # I think it is all managing advancing the state by a timestep because it's doing a lot of transport over bin boundaries
    # and I don't need to do any of that because of the continuous approach
    # growlg?
end
