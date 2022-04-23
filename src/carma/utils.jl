# 1.7 from Colella and Woodward
function advection_slope(dxjm1, dxj, dxjp1, ajm1, aj, ajp1)
    if xor(ajp1 < aj, aj < ajm1)
        return 0
    else
        prefactor = drj / (dxjm1 + drj + drjp1)
        diffterm = (2*dxjm1 + dxj) / (drjp1 + dxj) * (ajp1 - aj) + (dxj + 2dxjp1) / (dxjm1 + dxj) * (aj - ajm1)
        delta = prefactor * diffterm
        return min(abs(delta), abs(ajp1 - aj), abs(aj - ajm1)) * sign(delta)
    end
end

"""
1.6 from Colella and Woodward
convention of 4 elements per array, at indices [j-1, j, j+1, j+2]
"""
function advection_average(dx::AbstractVector, a::AbstractVector)
    dxjm1, dxj, dxjp1, dxjp2 = dx
    ajm1, aj, ajp1, ajp2 = a
    t1 = aj
    t2 = (dxj / (dxj + dxjp1)) * (ajp1 - aj)
    t3 = 1 / sum(dx)
    t4 = (2 * dxjp1 * dxj) / (dxj + dxjp1)
    t5 = (((dxjm1 + dxj) / (2 * dxj + dxjp1)) - (dxjp2 + dxjp1) / (2 * dxjp1 + dxj)) * (ajp1 - aj)
    t6 = -dxj * ((dxjm1 + dxj)/(2 * dxj + dxjp1))
    delp1 = advection_slope(dxj, dxjp1, dxjp2, aj, ajp1, ajp2)
    t7 = dxjp1 * ((dxjp1 + dxjp2) / (dxj + 2 * dxjp1))
    del = advection_slope(dxjm1, dxj, dxjp1, ajm1, aj, ajp1)
    return t1 + t2 + t3 * (t4 * t5 + t6 * delp1 + t7 * del)
end

function set_edges(aj_avg)
    
end

# move this to particle.jl once it's done
function dpdt_adv() 

end