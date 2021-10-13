using LinearAlgebra, BenchmarkTools, SuiteSparse

ACL = cholesky(mA)
@btime $x0 .= $ACL\$b

@btime CL = make_CL_in_julia(ACL);

x = similar(x0)
@btime ldiv_cl!(x,CL,b);

r1 = CartesianIndices(1:length(b))
r2 = CartesianIndices(CL.p)
xp = similar(x)
bp = similar(b)

@btime ldiv_cl!($x,$CL,$b,$r1,$r2,$xp, $bp);
@profiler @btime ldiv_cl!(x,CL,b,r1,r2,xp, bp);

sum(abs,x.-x0)


function ldiv_cl!(x,
                  CL::NamedTuple{(:L, :U, :p, :x_temp),
                      Tuple{SparseMatrixCSC{Float64, Int64},
                      SparseMatrixCSC{Float64, Int64},
                      Vector{Int64},
                      Matrix{Float64}}},
                  b, r1,r2,xp,bp)
    copyto!(bp,r1,b,r2);
    #@inbounds bp = view(b,CL.p)
    copyto!(xp,r1,x,r2);
    #@inbounds xp = view(x,CL.p)
    forward_substit!(CL.x_temp,CL.L, bp)
    backward_substit!(xp,CL.U, CL.x_temp)
    #invpermute!(x,CL.p)
end
