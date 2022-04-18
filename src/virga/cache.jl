"""
Cached variables that are virga-specific.
"""
struct VirgaCache
    fsed::Float64
    sig::Float64
    alpha_pressure::Pressure{Float64}

    function VirgaCache(
        fsed::Float64=0.5,
        sig::Float64=2.0,
        pressure::Vector{Pressure{Float64}},
        temperature::Vector{Temperature{Float64}},
        alpha_pressure::Union{Nothing,Pressure{Float64}}=nothing
    )
        if isnothing(alpha_pressure)
            alpha_pressure = minimum(pressure)
        end
        new(fsed, sig, alpha_pressure)
    end
end

