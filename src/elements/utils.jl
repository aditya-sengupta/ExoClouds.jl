# you may wonder why i bother having the subscripts if I'm just going to raise them, and to that I simply answer: The Aesthetic
"""
Raises the subscripts in element names.
"""
function unsubscribe(element::String)::String
    raiser = Dict("₂" => "2", "₃" => "3", "₄" => "4")
    map(x -> get(raiser, x, x), split(element, "")) |> join
end
