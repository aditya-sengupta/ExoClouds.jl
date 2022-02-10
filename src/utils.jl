using Unitful: 𝐌, 𝚯, 𝐋, 𝐓

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