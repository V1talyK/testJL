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
function repack(knLn, clL, rwL)
    fg = Vector(undef,length(knLn))

    for (k1, v1) in enumerate(knLn)
        fg[k1] = []
        for (k2,v2) in enumerate(v1)
            c1 = clL[v2]:clL[v2+1]-1
            push!(fg[k1],rwL[c1])
        end
    end

    fl = falses(length(knLn))
    for i = 2:length(knLn)
        fl[i] = length(intersect(vcat(fg[i]...),knLn[i-1]))>0
    end
    return fl
end

function test_kn(knLn)
    llvl = 1
    knLl = Int32.(vcat(knLn...))
    lvl_lng = Int32.(length.(knLn))
    iknL1 = Int32.(vcat(1,cumsum(length.(knLn)).+1)[1:end-1])
    iknL2 = Int32.(cumsum(length.(knLn)))
    knL32 = map(x->Int32.(x),knLn)

    knL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knL32[1])
    knL1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knLl)
    iknL1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=iknL1)
    lvl_lng_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=lvl_lng)

    sdf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.([length(iknL1),L.n, llvl]))
    cl.copy!(queue, y_bf, zr_bf);
    bs = (BLOCK_SIZE,)
    @time queue(krnL, bs, bs, zz_bf, knL1_bf, iknL1_bf, lvl_lng_bf,
                                y_bf, b_bf, clL_bf, rwL_bf, nzL_bf, sdf, lmem)
    yy = cl.read(queue, y_bf)

    e1r = [sum(abs,y0[knLn[i]].-yy[knLn[i]]) for i=1:length(knLn)]
        println(lineplot(e1r))
    println(findall(e1r.>100))
    return sum(abs,yy.-y0)
end
