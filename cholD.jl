CA = cholesky(-AS)
CA = SuiteSparse.CHOLMOD.cholesky(cuA)
SuiteSparse.CHOLMOD.lowrankupdate(CA, v1)

using SuiteSparse
