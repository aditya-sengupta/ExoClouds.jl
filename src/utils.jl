using Unitful: 𝐌, 𝚯, 𝐋, 𝐓

# Convention with questionable performance benefits (see if profiler throws a fit at this):
# dimensionless numbers are all type Real and are instantiated as whatever, most often Float64 but we'll let the compiler do its thing mostly
# (this is how Base Julia does it usually)
# quantities are all forced to be Float64s to allow for concrete allocations on a non-Base type

Temperature = Quantity{Float64,𝚯} # any number that carries units of temperature
Pressure = Quantity{Float64,𝐋𝐌𝐓⁻²}
Mass = Quantity{Float64, 𝐌}
Length = Quantity{Float64, 𝐋}
Acceleration = Quantity{Float64, 𝐋𝐓⁻²}

# you may wonder why i bother having the subscripts if I'm just going to raise them, and to that I simply answer: The Aesthetic
"""
Raises the subscripts in element names.
"""
function unsubscribe(element::String)::String
    raiser = Dict("₂" => "2", "₃" => "3", "₄" => "4")
    map(x -> get(raiser, x, x), split(element, "")) |> join
end