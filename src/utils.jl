using Unitful: 𝐌, 𝚯

Temperature = Quantity{<:Number,𝚯} # any number that carries units of temperature
Mass = Quantity{<:Number, 𝐌}