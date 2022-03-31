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
                      CL::NamedTuple{(:L, :U, :p, :x_temp, :y_temp),
                        Tuple{SparseMatrixCSC{Float32, Int64},
                        SparseMatrixCSC{Float32, Int64},
                        Vector{Int64},
                        Vector{Vector{Float32}},
                        Vector{Vector{Float32}}}},
                      b::Vector{Float32})
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


function forward_substit!(x, S, b)
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
        xc = xc/S.nzval[sc1]
        x[col] = xc
        setindex!(x,xc,col)
        sc1 = sc2
        #@inbounds sr=view(S.rowval,idx)

        for v in idx
             #foo2!(v,S,x,xc,sr,k)
             foo1!(v,S,x,xc)
        end
    end
end

function foo2!(v::Int64,
              S::SparseMatrixCSC{Float64, Int64},
              x::Vector{Float64},
              xc::Float64,
              sr,
              k)
    #y = S.nzval[v] * x[col];
    @inbounds z = S.nzval[v]
    @fastmath y = z * xc;
    #@inbounds i = S.rowval[v]
    @inbounds i = sr[k]
    #po = pointer(S.rowval)
    #i = unsafe_load(po, v)
    @inbounds z = getindex(x,i)
    z -= y
    @inbounds setindex!(x,z,i)
    return nothing
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
    @fastmath @inbounds for col = UU.n:-1:1
        idx = UU.colptr[col + 1]-2 :-1 : UU.colptr[col]
        #println(idx," ",S.rowval[idx])
        xc = x[col]
        xc = (b[col] + xc)/UU.nzval[UU.colptr[col+1]-1]
        x[col] = xc

        for v in idx
            foo1!(v,UU,x,xc)
        end
    end
end

function bar(v,UU,x,xc)
    x[UU.rowval[v]] -=  UU.nzval[v] * xc
end


y = zeros(Float32,CL.L.n)
function nme(y,S)
    @inbounds @fastmath @simd for col = 1:S.n
        y[col] = one(Float32)/S.nzval[S.colptr[col]]
    end
end

@btime nme($y,$CL32.L)
