
list = kn[2]
@time for (k,v) in enumerate(list)
    zz[k] = calc_zz_k(v,x,b,cl1,rw,nz)
end

s = css1(cl,rw,nz,x,row)

function calc_zz_k(row::Int64,
             x::Array{Float64,1},
             b::Array{Float64,1},
             cl::Array{Int64,1},
             rw::Array{Int64,1},
             nz::Array{Float64,1})

    sr1 = cl[row+1]-1

    s = css1(cl,rw,nz,x,row)
    #s1 = cssk(cl,rw,nz,x,row)
    #println(s-s1)
    @inbounds s = (b[row]-s)/nz[sr1];
    return s
end

function css1(cl::Array{Int64,1},rw::Array{Int64,1},nz::Array{Float64,1},
              x::Array{Float64,1},row)
    s = zero(Float64)

    sru = cl1[row]
    rng = sru:cl1[row+1]-2

    @inbounds @fastmath for v in rng
        z = getindex(rw,v)
        s += x[z]*nz[v]
    end
    #s = sdot1(x, nz, rw, rng)
    return s
end
