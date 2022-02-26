module Virga
    using Roots
    
    using Unitful: Mass, Pressure, Temperature
    using Unitful: bar, K
    using Interpolations

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

    function recommend_gas(pressure, temperature, mh, mmw)
        function choice(gas::Molecule)::Bool
            cond_p, cond_t = condensation_t(gas, mh, mmw, pressure)
            interp_cond_t = map(p -> fixed_interp(p, cond_p, cond_t), pressure)
            diff_curve = temperature .- interp_cond_t
            ((length(diff_curve[diff_curve .> 0*K]) > 0) && (length(diff_curve[diff_curve .< 0*K]) > 0))
        end
        filter(choice, available_molecules())
    end

    export condensation_t, recommend_gas
end # module