using Bijections

struct CARMACache end
    # I may really regret this design decision but we're going with it for now
    elements::Bijection{Int64,Element}

    function CARMACache(elements::AbstractArray{Element})
        new(Bijection(Dict(i => elements[i] for i in 1:length(elements))))
    end
end