function make_CL_in_julia(ACL, nth = 1)
    LL = sparse(ACL.L)
    UU = copy(LL')

    x_temp = [zeros(LL.n) for i in 1:nth]
    y_temp = [zeros(LL.n) for i in 1:nth]
    return (L = LL, U = UU, p = ACL.p, x_temp, y_temp)
end

function make_CL_in_julia32(ACL, nth = 1)
    LL = Float32.(sparse(ACL.L))
    UU = copy(LL')

    invS = zeros(Float32,LL.n)
    nme(invS,LL)

    x_temp = [zeros(Float32,LL.n) for i in 1:nth]
    y_temp = [zeros(Float32,LL.n) for i in 1:nth]
    return (L = LL, U = UU, p = ACL.p, x_temp, y_temp, invS)
end

function ldiv_cl!(x::Vector{Float64},
                      CL::NamedTuple{(:L, :U, :p, :x_temp, :y_temp),
                        Tuple{SparseMatrixCSC{Float64, Int64},
                        SparseMatrixCSC{Float64, Int64},
                        Vector{Int64},
                        Vector{Vector{Float64}},
                        Vector{Vector{Float64}}}},
                      b::Vector{Float64})
        @inbounds bp = view(b,CL.p)
        #@inbounds bp = b[CL.p]
        #@inbounds xp = view(x,CL.p)
        x_temp = CL.x_temp[Threads.threadid()]
        y_temp = CL.y_temp[Threads.threadid()]
        @inbounds copy!(y_temp,view(x,CL.p))

        forward_substit!(x_temp,CL.L, bp)
        backward_substit!(y_temp,CL.U, x_temp)
        @inbounds x[CL.p] .= y_temp
        return nothing
        #invpermute!(x,CL.p)
end

function ldiv_cl!(x::Vector{Float32},
                      CL::NamedTuple{(:L, :U, :p, :x_temp, :y_temp, :invS),
                        Tuple{SparseMatrixCSC{Float32, Int64},
                        SparseMatrixCSC{Float32, Int64},
                        Vector{Int64},
                        Vector{Vector{Float32}},
                        Vector{Vector{Float32}},
                        Vector{Float32}}},
                      b::Vector{Float32})
        @inbounds bp = view(b,CL.p)
        #@inbounds bp = b[CL.p]
        #@inbounds xp = view(x,CL.p)
        x_temp = CL.x_temp[Threads.threadid()]
        y_temp = CL.y_temp[Threads.threadid()]
        @inbounds copy!(y_temp,view(x,CL.p))

        forward_substit!(x_temp,CL.L, bp, CL.invS)
        backward_substit!(y_temp,CL.U, x_temp)
        @inbounds x[CL.p] .= y_temp
        return nothing
        #invpermute!(x,CL.p)
end


function forward_substit!(x, S, b, invS)
    x .= 0;
    sc1 = S.colptr[1]
    #sr = zeros(Int64,length(x))
    @inbounds for col = 1:S.n
        sc2 = S.colptr[col+1]
        idx = sc1+1 : sc2 - 1
        #idx = sc1+oneunit(sc1) : sc2 - oneunit(sc2)

        #println(idx," ",S.rowval[idx])
        xc = x[col]
        xc = (b[col] + xc)
        #xc = xc/S.nzval[sc1]
        xc = xc*invS[col]
        x[col] = xc
        setindex!(x,xc,col)
        sc1 = sc2
        #@inbounds sr=view(S.rowval,idx)

         for i in eachindex(idx)
             v = idx[i]
             z = foo3!(v,S,x,xc)
         end
    end
end

function foo2!(v::Int64,
              S::SparseMatrixCSC{Float64, Int64},
              x::Vector{Float64},
              xc::Float64)
    #y = S.nzval[v] * x[col];
    @inbounds z = S.nzval[v]
    @fastmath y = z * xc;
    @inbounds i = S.rowval[v]
    #@inbounds i = sr[k]
    #po = pointer(S.rowval)
    #i = unsafe_load(po, v)
    @inbounds z = getindex(x,i)
    @fastmath z -= y
    @inbounds setindex!(x,z,i)
    return nothing
end

function foo3!(v::Int64,
               S::SparseMatrixCSC{Float64, Int64},
               x::Vector{Float64},
               xc::Float64)
        #y = S.nzval[v] * x[col];
    @inbounds z = S.nzval[v]
    @fastmath y = z * xc;
    #@inbounds i = S.rowval[v]
    po = pointer(S.rowval)
    i = unsafe_load(po, v)
    @inbounds z = getindex(x,i)
    z -= y
    @inbounds setindex!(x,z,i)
    return z
end

function foo1!(v::Int64,
               S::SparseMatrixCSC{Float64, Int64},
               x::Vector{Float64},
               xc::Float64)
        #y = S.nzval[v] * x[col];
    @inbounds z = S.nzval[v]
    @fastmath y = z * xc;
    #@inbounds i = S.rowval[v]
    po = pointer(S.rowval)
    i = unsafe_load(po, v)
    @inbounds z = getindex(x,i)
    z -= y
    @inbounds setindex!(x,z,i)
    return nothing
end

function foo1!(v::Int64,
               S::SparseMatrixCSC{Float32, Int64},
               x::Vector{Float32},
               xc::Float32)
        #y = S.nzval[v] * x[col];
    @inbounds z = S.nzval[v]
    @fastmath y = z * xc;
    #@inbounds i = S.rowval[v]
    po = pointer(S.rowval)
    i = unsafe_load(po, v)
    @inbounds z = getindex(x,i)
    z -= y
    @inbounds setindex!(x,z,i)
    return nothing
end

function backward_substit!(x, UU, b)
    x .= 0;
    @inbounds @fastmath for col = UU.n:-1:1
        idx = UU.colptr[col + 1]-2 :-1 : UU.colptr[col]
        #println(idx," ",S.rowval[idx])
        xc = x[col]
        xc = (b[col] + xc)/UU.nzval[UU.colptr[col+1]-1]
        x[col] = xc
        for v in idx
            #foo1!(v,UU,x,xc)
            bar(v,UU,x,xc)
        end
    end
end

function bar(v,UU,x,xc)
    @inbounds @fastmath x[UU.rowval[v]] -=  UU.nzval[v] * xc
end


function nme(y,S)
    @inbounds @fastmath @simd for col = 1:S.n
        y[col] = one(Float32)/S.nzval[S.colptr[col]]
    end
end


function solv_krn1!(x::Array{Float64,1},
                   kn::Vector{Vector{Int64}},
                   b::Array{Float64,1},
                   zz::Array{Float64,1},
                   cl,rw,nz)

    list = kn[1]
    for (k,row) in enumerate(list)
        sr1 = cl[row+1]-1
        @inbounds zz[k] = b[row]/nz[sr1];
    end
    cp2!(x,zz,list)

    for i = 2:length(kn)
        list = kn[i]
        #calc_zz!(zz,kn[i],x,b,cl,rw,nz)
        #calc_zzt!(zz,kn[i],x,b,cl,rw,nz)
        calc_zz_kern(kernel,kn[i],zz,x,b,cl,rw,nz)
        cp2!(x,zz,list)
    end
end

function solv_krn!(x::Array{Float64,1},
                   kn::Vector{Vector{Int64}},
                   b::Array{Float64,1},
                   zz::Array{Float64,1},
                   cl,rw,nz)

    list = kn[1]
    for (k,row) in enumerate(list)
        sr1 = cl[row+1]-1
        @inbounds zz[k] = b[row]/nz[sr1];
    end
    cp2!(x,zz,list)

    for i = 2:length(kn)
        list = kn[i]
        calc_zz!(zz,kn[i],x,b,cl,rw,nz)
        #calc_zzt!(zz,kn[i],x,b,cl,rw,nz)
        #calc_zz_kern(kn[i],zz,x,b,cl,rw,nz)
        cp2!(x,zz,list)
    end
end

function calc_zz!(zz::Array{Float64,1},
                 list::Array{Int64,1},
                 x::Array{Float64,1},
                 b::Array{Float64,1},
                 cl::Array{Int64, 1},
                 rw::Array{Int64, 1},
                 nz::Array{Float64, 1})

    for (k,v) in enumerate(list)
        zz[k] = calc_zz_k(v,x,b,cl,rw,nz)
    end
end

function calc_zzt!(zz::Array{Float64,1},
                 list::Array{Int64,1},
                 x::Array{Float64,1},
                 b::Array{Float64,1},
                 cl::Array{Int64, 1},
                 rw::Array{Int64, 1},
                 nz::Array{Float64, 1})
    gh = collect(Iterators.partition(list,40))
    Threads.@threads for plist in gh
        for (k,v) in enumerate(plist)
            zz[k] = calc_zz_k(v,x,b,cl,rw,nz)
        end
    end
end

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

    sru = cl[row]
    rng = sru:cl[row+1]-2

    @inbounds @fastmath for v in rng
        z = getindex(rw,v)
        s += x[z]*nz[v]
    end
    #s = sdot1(x, nz, rw, rng)
    return s
end

function css(U::SparseMatrixCSC{Float64, Int64},x::Array{Float64,1},row::Int64)
    s::Float64 = 0.0
    sru = U.colptr[row]
    rng = sru:U.colptr[row+1]-2
    for v in rng
        z = getindex(U.rowval,v)
        @inbounds s += x[z]*U.nzval[v]
    end
    #s = sdot1(x, U.nzval, U.rowval, rng)
    return s
end


function make_order(U0)
    flag = true
    k = 1
    rn = []
    kn = []
    U = copy(U0)
    while flag
        rn = []
        for i=1:U.n
            if U.colptr[i]==U.colptr[i+1]-1
            #if length(L[i,:].nzind)==1
            #if count(i.==L.rowval)==1
                push!(rn,i)
            end
        end
        push!(kn,[])
        for i in rn
            push!(kn[k],i)
            U[i,:].=0
        end
        dropzeros!(U)
        k+=1
        flag = length(U.nzval)>0 & k<L.n
        println(k)
    end
    return map(x->Int64.(x), kn)
end
function cp2!(a,b,ind)
    for (k,v) in enumerate(ind)
        a[v] = b[k]
    end
end
