using SuiteSparse, LinearAlgebra

v = zeros(16610)
v[1]=1
lowrankdowndate(AS,v)
size(AS)

@time AS1 = SuiteSparse.CHOLMOD.lowrankupdate(AS,-v)
@time CL0=cholesky(-A);
@time for i=1:100 x0 = CL0\b; end;
lA = ldlt(-A)
@time for i=1:100 x1 = lA\b; end;

@time for i=1:100 lA1 = SuiteSparse.CHOLMOD.lowrankupdate(lA,v); end;
@time for i=1:100 CL1 = SuiteSparse.CHOLMOD.lowrankupdate(CL0,v); end;

sum(abs,x1.-x0)

x01 = lA1\b
x11 = CL1\b
sum(abs,x01-x11)
diag(lA1)[lA1.p.==1]

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

diag(lA).+diag(A)[lA.p]

y = lA.L\b[lA.p]
z = (sparse(1:16610,1:16610,diag(lA).+1))\y
x = lA.L'\z

x123 = (A-sparse(1:16610,1:16610,1))\b

x123[lA.p]
