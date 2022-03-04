using JLD, SparseArrays, BenchmarkTools, LinearAlgebra,CuthillMcKee
n = 1000
A = -sprand(n,n,0.001)
A = A + A'
rc = CartesianIndex.(1:n,1:n)
A[rc].=-sum(A,dims=2)[:].+1
p = symrcm(A)
A = A[:,p]
A = A[p,:]
b = zeros(n); b[1:4:end].=1
x0 = A\b


L1 = zeros(n,n)
n = A.n
L.=0
Li = zeros(n)
Lj = zeros(n)
a = zeros(n)
#Rdest = CartesianIndices(1:n)
println(1)
@time hcho(L,A,Li,Lj,n,a)
@btime hcho($L,$A,$Li,$Lj,$n,$a)
@profiler hcho(L,A,Li,Lj,n,a)
@code_warntype hcho(L,A,Li,Lj,n,a)

y = L'\b;
x1 = L\y;

sum(abs,x0.-x1)

@time cholesky(A);
@btime cholesky($A);

CI = CartesianIndices((1:10,1:10))
CI[:,1]

a = zeros(10)
b = rand(10,2)
Rdest = CartesianIndices((1:2:10,))
Rsrc = CartesianIndices((2:2:10,2:2))
copyto!(a,Rdest,vec(b),Rsrc)

A.colptr
A.rowval

r = zeros(Int64,0)
c = ones(Int64,1)
for i = 1:n
    v = A.colptr[i]:A.colptr[i+1]-1
    append!(r,minimum(A.rowval[v]):i)
    push!(c,length(r)+1)
end
#L = sprandn(n,n,0.1)
L = SparseMatrixCSC(n,n,c,r,zeros(length(r)))

VV = SparseVector(10,[1,3,7],[1,2,3])

L1 = L[:,j]
L1.


L.=0
@time hcho3(L,A,Li,Lj,n,a)

y = L'\b;
    x1 = L\y;
    sum(abs,x0.-x1)
