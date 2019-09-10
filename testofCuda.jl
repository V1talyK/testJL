n=16610;
elty = Float64
A1     = rand(elty,n,n); A[1,3] = 0;
A1     = sparse(A*A') #posdef
d_A   = CuArrays.CUSPARSE.CuSparseMatrixCSR(A);
b     = rand(elty,n)
d_b   = CuArray(b)
x     = zeros(elty,n)
d_x   = CuArray(x)
tol   = convert(real(elty),1e-4)
@time d_x   = CuArrays.CUSOLVER.csrlsvqr!(d_A,d_b,d_x,tol,one(Cint),'O')
h_x   = CuArrays.to_host(d_x)


A     = rand(elty,n,n)
A     = sparse(A*A') #posdef
d_A   = CuSparseMatrixCSR(-A);
b     = rand(elty,n)
d_b   = CuArray(b)
x     = zeros(elty,n)
d_x   = CuArray(x)
tol   = 10^2*eps(real(elty))
d_x   = CuArrays.CUSOLVER.csrlsvchol!(d_A,d_b,d_x,tol,zero(Cint),'O')

n=100;
A = sparse(rand(elty,n,n))
b = rand(elty,n)
x = zeros(elty,n)
tol = convert(real(elty),1e-6)
@time x =  CuArrays.CUSOLVER.csrlsvlu!(A,b,x,tol,one(Cint),'O')
@time d_x =  CuArrays.CUSOLVER.csrlsvlu!(d_A,d_b,d_x,tol,one(Cint),'O')
@time A\b;
