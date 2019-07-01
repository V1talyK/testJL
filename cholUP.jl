using SuiteSparse

v = ones(16610)
lowrankdowndate(AS,v)
size(AS)

@time AS1 = SuiteSparse.CHOLMOD.lowrankupdate(AS,-v)
@time CL0=cholesky(A);
@time x0 = CL0\b;

CL = copy(CL0)

v[2:end].=0;
v[1]=sqrt(2-A[1,1])
@time CL1 = SuiteSparse.CHOLMOD.lowrankupdate(CL,v)

x1 =CL1\b;

A2 = copy(A);
A2[1,1]=2
CL2=cholesky(A2);
x2 = CL2\b;

sum(abs.(x1-x2))
