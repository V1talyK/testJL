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
        for (k,row::Int64) in enumerate(list)
            # list[i]
            #ia = setdiff(L[row,:].nzind,row)
            #ia = U.rowval[U.colptr[row]:U.colptr[row+1]-2]
            #s = U.nzval[U.colptr[row]:U.colptr[row+1]-2]
            #x[row] = (b[row]-dot(view(x,ia),s))/L[row,row]

            sr1 = L.colptr[row]
            #s+= x[U.rowval[v]]*U.nzval[v]
            s = css(U,x,row)
            y = b[row]
            s = y-s
            @inbounds s = s/L.nzval[sr1]
            #@inbounds setindex!(x,s,row)
            #k+=1
            @inbounds zz[k] = s
        end
        x[list] .= view(zz,1:length(list))
    end
end

function css(U::SparseMatrixCSC{Float64, Int64},x::Array{Float64,1},row)
    s::Float64 = 0.0
    sru = U.colptr[row]
    rng = sru:U.colptr[row+1]-2
    @inbounds for v = rng
        z = getindex(U.rowval,v)
        s = s + x[z]*U.nzval[v]
    end
    return s
end


x2 = copy(x); x2.=0
@btime forward_substit!($x2, $L, $b, $invS)

invS = zeros(Float32,L.n)
nme(invS,L)

sum(abs,x0.-x2)
