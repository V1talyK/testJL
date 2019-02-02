using IterativeSolvers, Preconditioners, IncompleteLU

@time x0 = A\b;
@time x1 = IterativeSolvers.gmres(A, b);
@time x1 = IterativeSolvers.gmres(-A, b; Pl = p.L);
x11 = copy(x1)
@time IterativeSolvers.gmres!(x0, A, b; Pl = p.L);

@time x2 = IterativeSolvers.cg(-A, -b);
@time x2 = IterativeSolvers.cg(A, b; Pl = LUi);
@time x2 = IterativeSolvers.cg(A, b; Pl = LUf);
x21=zeros(Float64,length(x2));
x21=x0+20*rand(Float64,length(x2));
@time IterativeSolvers.cg!(x21, A, b; Pl = LU);



@time x3 = minres(A, b)
@time x3 = minres(A, b; Pl = p.L)

@time x4 = bicgstabl(A, b, 1; Pl = LUf);
@time x4 = bicgstabl(A, b, 1; Pl = LUi);
x44 = copy(x4)#+fill(1,length(x4[1]))
@time bicgstabl!(x44, A, b, 2 ; Pl = LU);

sum(abs.(x0-x1))
sum(abs.(x0-x2))
sum(abs.(x0-x4))

@time p = CholeskyPreconditioner(-A, 2)
@time p = AMGPreconditioner(A)
@time LUi = ilu(A, τ = 0.1)

x6 = @time jacobi(A, b; maxiter=5000)
sum(abs.(x0-x6))

@time x7 = gauss_seidel(A, b; maxiter=1000)

x8 = ssor(A, b, 0.01; maxiter=10000);
