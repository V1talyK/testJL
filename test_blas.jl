n=1000;
A = sprand(n,n,.05);
A = A+A';
A[1:n+1:n*n] = -sum(A,1);
A[1]=A[1]-1;
A[end]=A[end]-1;

B = rand(n);
B[1]=100;
B[end]=10;
@time x = (-A)\B;
@time x = A\B;

F = lufact(A);
L=F[:L]
U=F[:U]
p = F[:p]
y = BLAS.trsv('l', 'n', 'u', L, B)

Bp=copy(B)
Bp[p]=B;
y=L\B[p];
x1 =U\y;
x2=copy(x1);
x2[p]=x1;


@time CH = cholfact(-A);
@time x1 = CH\B;
