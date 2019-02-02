using IterativeSolvers, Preconditioners

@time x0 = A\b;
@time x1 = IterativeSolvers.gmres(A, b);
@time x11 = IterativeSolvers.gmres!(x0, A, b);

@time x2 = IterativeSolvers.cg(A, b);
x21=zeros(Float64,length(x2));
x21=x0+1#rand(Float64,length(x2));
@time IterativeSolvers.cg!(x21, A, b);


p = CholeskyPreconditioner(A, 2)

@time x3 = minres(A, b)

sum(abs.(x0-x11))
sum(abs.(x0-x2))
