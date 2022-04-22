function setup(nbin::Int64, nelem::Int64, ngas::Int64)
    state_init = State(zeros(nbin,nelem), zeros(ngas), 0.0K)
    
    return state_init
end