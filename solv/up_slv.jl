x0 = ACL\b
y = zeros(length(b))
bp = b[ACL.p]
y0 = L\bp
@time solv_krn!(y,kn,bp,zz, cl1,rw,nz)

L0=sparse(ACL.L)
U0 = copy(sparse(L0)')

xp = similar(x0)

knU = make_order1(L0)
clU = copy(L0.colptr)
rwU = copy(L0.rowval)
nzU = copy(L0.nzval)

knL = make_order(U0)
clL = copy(U0.colptr)
rwL = copy(U0.rowval)
nzL = copy(U0.nzval)

xp = zeros(size(x0))
@time solv_krnU!(xp,kn1,y,zz, clL,rwL,nzL)
x[ACL.p] = xp

@time like_gpu_slv!(x, bp, zz, clU, rwU, nzU, clL, rwL, nzL, knL, knU)

sum(abs,x.-x0)
sum(abs,y.-y0)

@time ACL\b;

function like_gpu_slv!(x, bp, zz, clU, rwU, nzU, clL, rwL, nzL, knL, knU)
    solv_krn!(y,knL,bp,zz, clL,rwL,nzL)
    solv_krnU!(xp,knU,y,zz, clU,rwU,nzU)
    x[ACL.p] .= xp;
    return nothing
end


function solv_krnU!(x::Array{Float64,1},
                   kn::Vector{Vector{Int64}},
                   b::Array{Float64,1},
                   zz::Array{Float64,1},
                   cl,rw,nz)

    list = kn[1]
    for (k,row) in enumerate(list)
        sr1 = clL[row+1]-1
        @inbounds zz[k] = b[row]/nz[sr1];
    end
    cp2!(x,zz,list)

    for i = 2:length(kn)
        list = kn[i]
        calc_zzU!(zz,kn[i],x,b,clL,rwL,nzL)
        #calc_zzt!(zz,kn[i],x,b,cl,rw,nz)
        #calc_zz_kern(kn[i],zz,x,b,cl,rw,nz)
        cp2!(x,zz,list)
    end
end

function calc_zzU!(zz,
                 list,
                 x::Array{Float64,1},
                 b::Array{Float64,1},
                 cl::Array{Int64, 1},
                 rw::Array{Int64, 1},
                 nz::Array{Float64, 1})

    for (k,v) in enumerate(list)
        zz[k] = calc_zz_kU(v,x,b,cl,rw,nz)
        #put!(mcnl3,(v,k))
    end
end

function calc_zz_kU(row::Int64,
             x::Array{Float64,1},
             b::Array{Float64,1},
             cl::Array{Int64,1},
             rw::Array{Int64,1},
             nz::Array{Float64,1})

    sr1 = cl[row]

    s = cssU(cl,rw,nz,x,row)
    #s1 = cssk(cl,rw,nz,x,row)
    #println(s-s1)
    @inbounds s = (b[row]-s)/nz[sr1];
    return s
end

function cssU(cl::Array{Int64,1},rw::Array{Int64,1},nz::Array{Float64,1},
              x::Array{Float64,1},row)
    s = zero(Float64)

    sru = cl[row]+1
    rng = sru:cl[row+1]-1

    @inbounds @fastmath for v in rng
        z = getindex(rw,v)
        s += x[z]*nz[v]
    end
    #s = sdot1(x, nz, rw, rng)
    return s
end
