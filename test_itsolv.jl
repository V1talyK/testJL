using IterativeSolvers, Preconditioners, IncompleteLU, LinearAlgebra, BenchmarkTools
using SuiteSparse, SparseArrays
LinearAlgebra.BLAS.set_num_threads(1) #set_num_threads

@btime x0 = A\b;
@btime ACL = cholesky(mA)
@btime x0 = ACL\b;

@btime ALU = lu(mA)
@btime x0 = ALU\b;


@btime x2 = IterativeSolvers.cg(mA, b);
x2.=0
@btime IterativeSolvers.cg!($x2, $mA, $b);
@btime IterativeSolvers.cg!(x2, mA, b; maxiter = 1, log = true, Pl = LUi);
@btime IterativeSolvers.cg!(x2, mA, b; maxiter = 100, log = true);
@btime IterativeSolvers.cg!($x2, $mA, $b; Pl = $LUi);

@time for i=1:50
    x0lu = LUf\b;
end

@time x1 = IterativeSolvers.gmres(A, b);
@time x1 = IterativeSolvers.gmres(-A, b; Pl = p.L);
x11 = copy(x1)
@time IterativeSolvers.gmres!(x0, A, b; Pl = p.L);


@btime IterativeSolvers.jacobi!(x2, mA, b; maxiter=10)



@time x2 = IterativeSolvers.cg(A, b; Pl = LUi);

@btime IterativeSolvers.cg!($x2,$mA, $b; Pl = $LUi);
@profiler IterativeSolvers.cg!(x2,mA, b; Pl = LUi);
IterativeSolvers.cg!(x2,mA, b; Pl = LL);

@time for i=1:50s
    x3 = IterativeSolvers.cg(A, b; Pl = LUi);
end

cux2=CuArrays.CuArray(0*rand(length(x)))

@time IterativeSolvers.cg!(cux2,cuA, cuB);
@time cux3 = IterativeSolvers.cg(cuA, cuB)

cuL =  CuArrays.CUSPARSE.CuSparseMatrixCSC(LUi.L);
cux3 = IterativeSolvers.cg(cuA, cuB; Pl = cuL);

@time for i=1:50
    #x2 = IterativeSolvers.cg(A, b; Pl = LUi);
    #x2 = IterativeSolvers.cg(A, b);
    cux3 = IterativeSolvers.cg(cuA, cuB; Pl = iluA.L);
end

@time for i=1:50
    x2 = IterativeSolvers.cg(A, b; Pl = LUf);
end

@time x2 = IterativeSolvers.cg(A, b; Pl = LUf);
x21=zeros(Float64,length(x2));
x21=x0+20*rand(Float64,length(x2));
@time IterativeSolvers.cg!(x21, A, b; Pl = LU);

@time cuA*d_b;

p = CholeskyPreconditioner(-cuA, 2)

@time x3 = minres(A, b)
@time x3 = minres(A, b; Pl = LUi.L)

@time x4 = bicgstabl(A, b, 1; Pl = LUi);
@time x4 = bicgstabl(cuA, cuB, 1);
@time for i=1:50
    x4 = bicgstabl(A, b, 1; Pl = LUi);
end
x44 = copy(x4)#+fill(1,length(x4[1]))
@time bicgstabl!(x44, A, b, 1 ; Pl = LUi);

sum(abs.(x0-x11))
sum(abs.(x0-x2))
sum(abs.(x0-x4))

@time p = CholeskyPreconditioner(mA, 2)
@time p = AMGPreconditioner(A)
@time LUi = ilu(mA, τ = 0.000001)

x6 = @time jacobi(A, b; maxiter=1000, Pl = LUi)
sum(abs.(x0-x6))

@time x7 = gauss_seidel(A, b; maxiter=1000)

x8 = ssor(A, b, 0.01; maxiter=10000);
