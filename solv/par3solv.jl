using SparseArrays, BenchmarkTools, LinearAlgebra, LoopVectorization

r = vcat(collect(1:9),4:9,8:9)
c = vcat(1:9,1,1:5,5:-1:4)
L = sparse(r,c,1)
f = collect(0:8)
x0 = L\f


A = sprand(100,100,0.1)
A = A+A';
A[CartesianIndex.(1:100,1:100)].=-sum(A,dims=1)[:].-1
A = -A
ACL = cholesky(mA)
L = sparse(ACL.L)
b = rand(100)
x0 = L\b

flag = true
k = 1
rn = []
kn = []
while flag
    rn = []
    for i=1:L.n
        if length(L[i,:].nzind)==1
        #if count(i.==L.rowval)==1
            push!(rn,i)
        end
    end
    push!(kn,[])
    for i in rn
        push!(kn[k],i)
        L[:,i].=0
    end
    dropzeros!(L)
    k+=1
    flag = length(L.nzval)>0 & k<L.n
    println(k)
end

U = copy(L')
x = zeros(length(b))
@time solv_krn!(x,kn1,b,L,U)
@btime solv_krn!($x,$kn1,$b,$L,$U)
@btime $x0.=$L\$b
@profiler for i=1:500 solv_krn!(x,kn1,b,L,U); end;

sum(abs,x.-x0)
kn1 = map(x->Int64.(x), kn)

function solv_krn!(x::Array{Float64,1},kn::Vector{Vector{Int64}},b::Array{Float64,1},L::SparseMatrixCSC{Float64, Int64},U::SparseMatrixCSC{Float64, Int64})
    zz = zeros(length(b))
    for j = 1:length(kn)
        list = kn[j]
        #k=0
        @inbounds zz[1:length(list)] .= 1.0 ./ view(L.nzval,view(L.colptr,list))
        for (k,row::Int64) in enumerate(list)
            dfg!(zz,L,U,x,b,row)
        end
        x[list] .= view(zz,1:length(list))
    end
end

function dfg!(zz,L,U,x,b,row,k)

    # list[i]
    #ia = setdiff(L[row,:].nzind,row)
    #ia = U.rowval[U.colptr[row]:U.colptr[row+1]-2]
    #s = U.nzval[U.colptr[row]:U.colptr[row+1]-2]
    #x[row] = (b[row]-dot(view(x,ia),s))/L[row,row]

    sr1 = L.colptr[row]
    sru = U.colptr[row]
    sru2 = U.colptr[row+1]-2
    rng = sru:sru2
    s = css(rng,U,x)

    y = b[row]
    s = y-s
    @inbounds @fastmath s = s*zz[k]
    #@inbounds setindex!(x,s,row)
    #k+=1
    @inbounds zz[k] = s
end

function css(rng, U::SparseMatrixCSC{Float64, Int64},x::Array{Float64,1})
    s::Float64 = 0.0
    for v in rng
        #s+= x[U.rowval[v]]*U.nzval[v]
        @inbounds z = getindex(U.rowval,v)
        @inbounds s = s + x[z]*U.nzval[v]
    end
    return s
end


x2 = copy(x); x2.=0
@time forward_substit!(x2, L, b, invS)
@btime forward_substit!($x2, $L, $b, $invS)

invS = zeros(Float32,L.n)
nme(invS,L)

sum(abs,x0.-x2)
