using RecursiveFactorization
using IncompleteLU
using ILUZero
using LimitedLDLFactorizations

LU = RecursiveFactorization.lu(A); # out-of-place

SparseArrays.UMFPACKFactorization(;reuse_symbolic=true)
SuiteSparse.lufact(mA; reuse_symbolic=true)

AILU = ilu(A, τ = 0.001)
xilu = AILU\b

AILU = ilu(A, τ = 0.001)
xilu = AILU\b

LLDL = lldl(A, memory = 5)
xldl = LLDL\b


@btime IterativeSolvers.cg!($xilu, $A, $b; Pl = AILU, reltol=1e-2);
IterativeSolvers.cg!(xilu, A, b; Pl = AILU, reltol=1e-2);

mean(abs, x0.+xilu)
extrema(x0.+xilu)

@btime IterativeSolvers.cg!($xldl, $A, $b; Pl = LLDL, reltol=1e-2);
IterativeSolvers.cg!(xldl, A, b; Pl = LLDL, reltol=1e-2);

mean(abs, x0.+xldl)
extrema(x0.+xldl)
