using IterativeSolvers, Preconditioners

@time x0 = A\b;
@time x1 = IterativeSolvers.gmres(A, b);
@time x1 = IterativeSolvers.gmres(A, b; Pl = p.L);
x11 = copy(x1)
@time IterativeSolvers.gmres!(x0, A, b; Pl = p.L);

@time x2 = IterativeSolvers.cg(A, b; Pl = p.L);
x21=zeros(Float64,length(x2));
x21=x0+1#rand(Float64,length(x2));
@time IterativeSolvers.cg!(x21, A, b);



@time x3 = minres(A, b)

sum(abs.(x0-x1))
sum(abs.(x0-x2))

@time p = CholeskyPreconditioner(-A, 2)
