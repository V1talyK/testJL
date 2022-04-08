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
zz = zeros(maximum(length.(kn)))

x = zeros(length(b))
@time solv_krn!(x,kn,b,L,U,zz, fg)
@btime solv_krn!($x,$kn,$b,$L,$U,$zz,$fg)
@btime $x0.=$L\$b
@profiler for i=1:500 solv_krn!(x,kn,b,L,U,zz, fg); end;
@profile solv_krn!(x,kn,b,L,U, zz, fg)

sum(abs,x.-x0)

cl = copy(U.colptr)
rw = copy(U.rowval)
nz = copy(U.nzval)
fg(y::Int64) = fgh(y::Int64,L,U,x,b,cl,rw,nz)

list = kn[1]
@profiler bbr(zz,fg,list)

function bbr(zz::Array{Float64,1},fg::Function,list::Array{Int64,1})
    for (k,v) in enumerate(list)
        zz[k] = fg(v)
    end
end

x2 = copy(x); x2.=0
@time forward_substit!(x2, L, b, invS)
@btime forward_substit!($x2, $L, $b, $invS)

invS = zeros(Float32,L.n)
nme(invS,L)

sum(abs,x0.-x2)
@code_warntype solv_krn!(x,kn1,b,L,U)
@code_warntype fg(kn[1][1])


y = zeros(Float32,CL.L.n)
@btime nme($y,$CL32.L)
