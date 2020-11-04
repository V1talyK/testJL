using Pkg, SparseArrays
ENV["JULIA_PARDISO"]="/home/lik/temp"
ENV["MKLROOT"] = "/opt/intel/compilers_and_libraries_2019.0.117/linux/mkl"
ENV["PARDISO_LIC_PATH"] = "/home/lik/temp"
Pkg.build("Pardiso")
using Pardiso
Pardiso.show_build_log()

Pkg.add("MKLSparse")

using MKL, MKLSparse, LinearAlgebra
BLAS.vendor()

ps = PardisoSolver()

A = sparse(rand(10, 10))
B = rand(10)
X = zeros(length(b))
Pardiso.solve(ps, A, B)
solve(ps, A, B)
@time A\B



verbose = false

m = 3  # The number of right-hand sides.
n = 4  # The number of equations.

A = sparse([ 0. -2  3 0
            -2  4 -4 1
            -3  5  1 1
             1 -3  0 2 ])

# Generate a random collection of right-hand sides.
B = ones(n,m)

# Initialize the PARDISO internal data structures.
ps = PardisoSolver()

if verbose
    set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
end

# If we want, we could just solve the system right now.
# Pardiso.jl will automatically detect the correct matrix type,
# solve the system and free the data
X1 = solve(ps, A, B)
