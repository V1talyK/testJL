using Krylov, ILUZero, CUDA, LimitedLDLFactorizations
@btime xK, stats = Krylov.cg(A, b)
@btime x = A\b;


perm = 1:A.n
@btime LLDL = lldl(mA, perm, memory = 5)
B = tril(A)
LLDL = lldl(B, perm, memory = 5)
@btime
y = LLDL\b
x = LLDL.L\y


@btime LU = ilu0(mA)

ldiv!(x, LU, b)
@btime ldiv!($x, $LU, $b)

x2 = similar(x)

@btime x2 = IterativeSolvers.cg(mA, b);
x2.=0
@btime IterativeSolvers.cg!($x2, $mA, $b);
@btime IterativeSolvers.cg!($x2, $mA, $b; Pl = $LU);




A_cpu = rand(20, 20)
b_cpu = rand(20)

# GPU Arrays
A_gpu = CuMatrix(A_cpu)
b_gpu = CuVector(b_cpu)

# Solve a square and dense system on GPU
x, stats = bilq(A_gpu, b_gpu)

A_cpu = sprand(200, 200, 0.3)
b_cpu = rand(200)

# GPU Arrays
A_gpu = CUDA.CUSPARSE.CuSparseMatrixCSC(A_cpu)
b_gpu = CuVector(b_cpu)

# Solve a rectangular and sparse system on GPU
x, stats = Krylov.lsmr(A_gpu, b_gpu)
A_gpu*b_gpu
