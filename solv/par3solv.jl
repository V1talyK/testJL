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
U = copy(L')
kn = make_order(U)

x = zeros(length(b))
@time solv_krn!(x,kn,b,L,U,zz)
@btime solv_krn!($x,$kn,$b,$L,$U)
@btime $x0.=$L\$b
@profiler for i=1:500 solv_krn!(x,kn,b,L,U,zz); end;
@profile solv_krn!(x,kn1,b,L,U)

sum(abs,x.-x0)
zz = zeros(maximum(length.(kn)))

fg(y) = fgh(y,L,U,x,b)
function solv_krn!(x::Array{Float64,1},kn::Vector{Vector{Int64}},b::Array{Float64,1},L::SparseMatrixCSC{Float64, Int64},U::SparseMatrixCSC{Float64, Int64},zz::Array{Float64,1})

    for j = 1:length(kn)
        list::Vector{Int64} = kn[j]
        for (k,v) in enumerate(list)
            s = fg(v)
            zz[k] = s
        end
    end
end

function css(U::SparseMatrixCSC{Float64, Int64},x::Array{Float64,1},row::Int64)
    s::Float64 = 0.0
    sru = U.colptr[row]
    rng = sru:U.colptr[row+1]-2
    for v in rng
        z = getindex(U.rowval,v)
        @inbounds s += x[z]*U.nzval[v]
    end
    return s
end


x2 = copy(x); x2.=0
@time forward_substit!(x2, L, b, invS)
@btime forward_substit!($x2, $L, $b, $invS)

invS = zeros(Float32,L.n)
nme(invS,L)

sum(abs,x0.-x2)
@code_warntype solv_krn!(x,kn1,b,L,U)


function fgh(row,L::SparseMatrixCSC{Float64, Int64},U::SparseMatrixCSC{Float64, Int64},x::Array{Float64,1},b::Array{Float64,1})
    #row = list[i]
    #ia = setdiff(L[row,:].nzind,row)
    #ia = U.rowval[U.colptr[row]:U.colptr[row+1]-2]
    #s = U.nzval[U.colptr[row]:U.colptr[row+1]-2]
    #x[row] = (b[row]-dot(view(x,ia),s))/L[row,row]

    sr1 = L.colptr[row]
    #s+= x[U.rowval[v]]*U.nzval[v]
    sru = U.colptr[row]
    rng = sru:U.colptr[row+1]-2
    s::Float64 = css(U,x,row)
    #s = dot(view(x,view(U.rowval,rng)),view(U.nzval,rng))
    #s =
    @inbounds s = (b[row]-s)/L.nzval[sr1];
    #@inbounds x[row] = s
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
