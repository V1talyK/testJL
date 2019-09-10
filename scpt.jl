#ENV["MUMPS_PREFIX"] = "/usr/lib/x86_64-linux-gnu/"
#ENV["SCALAPACK_PREFIX"] = "/usr/lib/x86_64-linux-gnu/"
ENV["OMP_NUM_THREADS"]=4
using MUMPS
using MPI
using SparseArrays, LinearAlgebra, JLD
d = JLD.load(joinpath("/home/lik/Документы/proto/testJL","myfile200k.jld"));
A=d["A"];
b=d["b"];

MPI.Init()
mumps = Mumps{Float64}(mumps_unsymmetric, default_icntl, default_cntl64)  # Real, general unsymmetric

et=@time begin
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
println(ENV["OMP_NUM_THREADS"])
#println(et)
#
