n = 5000;
A = rand(n,n);
B = rand(n);

cA = cu(A);
cB= cu(B);

@time C = A*B;
@time cC = cA*cB;

@time for i=1:10
    C = A*B;
end

@time for i=1:10
    cC = cA*cB;
end

@test C ≈ Array(cC)

As = sprand(n,n,0.005);
As+=As';
A1=As-sparse(1:n,1:n,sum(As,dims = 2)[:])

cA1 = cu(A1)
L_U = lu(A1);
cL_U = lu(cA1);

d_A = CuArrays.CUSOLVER.potrf!('U',cA)

A    = A*A';
A = SparseArrays.sparse(A);
d_A  = CuArrays.CuArray(A1)
d_A1  = CuArrays.CUSOLVER.potrf!('U',d_A)

d_A = CuArrays.CuArray(A1)
d_A1 = similar(d_A);
d_B = CuArrays.CuArray(B)
d_C = similar(d_B);
@time begin
    d_A1,d_ipiv = CuArrays.CUSOLVER.getrf!(d_A)
    d_C = CuArrays.CUSOLVER.getrs!('N',d_A1,d_ipiv,d_B)
end
h_B = collect(d_B)

sum(h_B - C1)
@time C1 = A1\B

CuArrays.CUSOLVER.csrlsvlu(A1)
CuArrays.CUSOLVER.op(A1)
