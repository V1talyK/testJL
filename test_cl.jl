using CLBlast, OpenCL
@static if VERSION < v"0.7-"
    LA = LinAlg
else
    using Random, LinearAlgebra
    LA = LinearAlgebra
end

device, context, queue = cl.create_compute_context()

# setup data
α = 1.f0
β = 1.f0
n=1
A = rand(Float32, 100*n, 100*n)
B = rand(Float32, 100*n, 1)
C = zeros(Float32, 100*n, 60*n)

# transfer data
A_cl = cl.CLArray(queue, A)
B_cl = cl.CLArray(queue, B[:])
C_cl = cl.CLArray(queue, C)

# compute
@time LA.BLAS.gemm!('N', 'N', α, A, B, β, C)
@time CLBlast.gemm!('N', 'N', α, A_cl, B_cl, β, C_cl)

# compare results
Array(C_cl)
@assert cl.to_host(C_cl) ≈ C


@time fA\b;
@time A\b;

@time A*B;
@time A_cl*B_cl;
