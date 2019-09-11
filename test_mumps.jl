ENV["MUMPS_PREFIX"] = "/usr/lib/x86_64-linux-gnu/"
ENV["OMP_NUM_THREADS"]=4
using MUMPS
using MPI
using SparseArrays, LinearAlgebra
MPI.Init()
A = sprand(10, 10, .2)+sparse(1:10,1:10,1); rhs = rand(10)
x =MUMPS.solve(mumps, rhs)  # Mumps object is created and destroyed
norm(x - A \ rhs) / norm(x)
MPI.Finalize()     # if you'

mumps = Mumps{Float64}(mumps_unsymmetric, default_icntl, default_cntl64)  # Real, general unsymmetric
A = sparse(rand(4,4)); rhs = rand(4)       # Happens on all cores

b = copy(rhs)
@time begin
    for i=1:2
        associate_matrix!(mumps, A)
        factorize!(mumps)
        associate_rhs!(mumps, b)
        solve!(mumps)
        x = get_solution(mumps)
    end
end
finalize(mumps)
MPI.Finalize()


@time A\b;
