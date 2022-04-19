function condensation_t(
    atm::Atmosphere,
    gas::Element, 
    pressures::Union{Nothing,Vector{Pressure{Float64}}}=nothing
)
    if isnothing(pressures)
        pressures = each(atm.P)
    end
    pressures, [
        find_zero(
                t -> find_cond_t(t, atm.P, atm.mh, atm.mmw, gas),
                (10.0*K, 10000.0*K),
                Roots.Brent()
            )
        for p in pressures
    ]
end

function recommend_gas(atm::Atmosphere, temperature::Vector{Temperature{Float64}}; makeplot=true)
    p = plot(yaxis=:log10, yflip=true, xlabel="Temperature (K)", ylabel="Pressure (bar)", yticks=10.0 .^(-5:2))
    if makeplot
        plot!(temperature ./ K, pressure ./ bar, label="User", linestyle=:dash, size=(800,600))
    end

    function choice(gas::Element)::Bool
        cond_p, cond_t = condensation_t(atm, gas)
        interp_cond_t = (LinearInterpolation((cond_p,), cond_t, extrapolation_bc=Flat())).(each(atm.P))
        diff_curve = temperature .- interp_cond_t
        rec = (maximum(diff_curve) > 0.0K) && (minimum(diff_curve) < 0.0K)
        if makeplot
            width = rec ? 5 : 1
            plot!(cond_t ./ K, cond_p ./ bar, lw=width, label=typeof(gas))
        end
        rec
    end

    filter(choice, available_elements()), p
end