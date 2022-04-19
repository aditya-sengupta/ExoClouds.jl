"""
Cached variables that are virga-specific.

Note: the pressure, temperature etc that get passed into the constructor here are not necessarily the same as the ones passed into Atmosphere. They can be, but we aren't going to structurally assume that, because atmospheres are continuous and simulations run on top of them are discrete. So we'll primarily interface with pressure, temperature etc. as functions P(z), T(z), etc.

I don't like that this is so grid dependent — I'd prefer to have continuous functions that are evaluated on the fly. But for now I just need something that works here.
"""
struct VirgaCache
    fsed::Float64
    sig::Float64
    alpha_pressure::Pressure{Float64}
    temperature::Vector{Temperature{Float64}}
    layer_temperature::Vector{Temperature{Float64}}
    layer_pressure::Vector{Pressure{Float64}}
    kz::Vector{Float64}
    layer_dz::Vector{Length{Float64}}
    z_alpha::Length{Float64}

    function VirgaCache(
        atm::Atmosphere,
        fsed::Float64=0.5,
        sig::Float64=2.0,
        pressure::Vector{Pressure{Float64}},
        temperature::Vector{Temperature{Float64}},
        alpha_pressure::Union{Nothing,Pressure{Float64}}=nothing,
        kz::Union{Nothing,Vector{Float64}}=nothing,
        chf::Union{Nothing,Vector{Float64}}=nothing;
        convective_overshoot = 1/3
    )
        if isnothing(alpha_pressure)
            alpha_pressure = minimum(pressure)
        end
        dlnp = log.(pressure[2:end] ./ pressure[1:end-1])
        dtdlnp = diff(temperature) / dlnp
        layer_pressure = (1/2) .* (pressure[2:end] .+ pressure[1:end-1])
        layer_temperature = temperature .+ log(pressure[2:end] ./ layer_pressure) .* dtdlnp
        H = scale_height.(atm, layer_temperature)
        layer_dz = H .* dlnp
        mixl = mixing_length.(atm, layer_temperature, layer_pressure)
        if isnothing(kz)
            @assert !isnothing(chf) "need to specify chf or kz"
            for iz = (length(layer_pressure):-1:0)
                ratio_min = convective_overshoot .* pressure[iz] / pressure[iz+1]
                chf[iz] = max(chf[iz], ratio_min * chf[iz+1])
            end
            # vertical eddy diffusion coefficient (cm^2/s)
            # from Gierasch and Conrath (1985)
            r = R / atm.mmw
            gc_kzz = (1/3) .* H .* (mixl ./ H) .^ (4/3) * (r * chf[2:end] / (atm.rho(atm.zref) * atm.cₚ)) .^ (1/3)
            kz = maximum.(gc_kzz, kz_min)            
        else
            kz = 0.5 * (kz[2:end] .+ kz[1:end-1])
        end

        p_alpha = findnearest(layer_pressure, alpha_pressure)
        z_temp = layer_dz |> reverse |> cumsum |> reverse
        z_alpha = z_temp[p_alpha]
        
        new(fsed, sig, alpha_pressure, temperature, layer_temperature, layer_pressure, kz, layer_dz, z_alpha)
    end
end

# get_kz_mixl

p_alpha = find_nearest_1d(self.p_layer/1e6, self.alpha_pressure)
z_temp = np.cumsum(self.dz_layer[::-1])[::-1]
self.z_alpha = z_temp[p_alpha]