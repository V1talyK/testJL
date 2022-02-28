using JLD, SparseArrays, BenchmarkTools, LinearAlgebra
n = 1000
A = -sprand(n,n,0.2)
A = A + A'
rc = CartesianIndex.(1:n,1:n)
A[rc].=-sum(A,dims=2)[:].+1

b = zeros(n); b[1:4:end].=1
x0 = A\b


L = zeros(n,n)
Li = zeros(n)
Lj = zeros(n)
a = zeros(n)
#Rdest = CartesianIndices(1:n)

@time hcho(L,A,Li,Lj,n,a)
@btime hcho($L,$A,$Li,$Lj,$n,$a)
@profiler hcho(L,A,Li,Lj,n,a)

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
