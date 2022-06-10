function get_max_wide(knLn,clL)
    sd = Vector(undef,length(knLn))
    for (k,kni) in enumerate(knLn)
        sd[k] = []
        for trow in kni
            push!(sd[k],length(clL[trow]:clL[trow+1]-1))
        end
    end

    ia = argmax.(sd)
    ib = argmax(getindex.(sd,ia))
    return maximum(maximum.(sd)), (ib, ia[ib])
end

function cut_lvl_by_BS(knL,BLOCK_SIZE)
    knLn = Vector(undef,0)
    lvl = Vector(undef,0)
    for (k,v) in enumerate(knL)
        for (k1,v1) = enumerate(Iterators.partition(v,BLOCK_SIZE))
            push!(knLn,v1)
            push!(lvl,k)
        end
    end
    return knLn, lvl
end
