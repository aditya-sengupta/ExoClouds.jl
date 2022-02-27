module Virga
    using Interpolations
    using Roots
    using Plots

    using Unitful: Mass, Pressure, Temperature
    using Unitful: bar, K

    using ..Molecules

    include("roots.jl")
    include("utils.jl")
    export find_cond_t

    function condensation_t(
        gas::Molecule, 
        mh::Real, 
        mmw::Mass, 
        pressure::Vector{<:Pressure}=(10 .^(-6:8/19:2)) * bar
    )
        pressure, Vector{Temperature}([
            find_zero(
                    t -> find_cond_t(t, p, mh, mmw, gas),
                    (10*K, 10000*K),
                    Roots.Brent()
                )
            for p in pressure
        ])
    end

    function recommend_gas(pressure::Vector{<:Pressure}, temperature::Vector{<:Temperature}, mh::Real, mmw::Mass; makeplot=true)
        p = plot(yaxis=:log10, yflip=true, xlabel="Temperature (K)", ylabel="Pressure (bar)", yticks=10.0 .^(-5:2))
        if makeplot
            plot!(temperature ./ K, pressure ./ bar, label="User", linestyle=:dash, size=(800,600))
        end

        function choice(gas::Molecule)::Bool
            cond_p, cond_t = condensation_t(gas, mh, mmw)
            interp_cond_t = map(p -> fixed_interp(p, cond_p, cond_t), pressure)
            diff_curve = temperature .- interp_cond_t
            rec = (maximum(diff_curve) > 0K) && (minimum(diff_curve) < 0K)
            if makeplot
                width = rec ? 5 : 1
                plot!(cond_t ./ K, cond_p ./ bar, lw=width, label=typeof(gas))
            end
            rec
        end

        filter(choice, available_molecules()), p
    end

    export condensation_t, recommend_gas
end # module