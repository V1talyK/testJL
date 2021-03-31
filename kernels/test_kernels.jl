using KernelAbstractions, Test, CUDA, LinearAlgebra

if has_cuda_gpu()
    CUDA.allowscalar(false)
end

# Simple kernel for matrix multiplication
@kernel function matmul_kernel!(a, b, c)
    i, j = @index(Global, NTuple)
    #println(i,j)
    # creating a temporary sum variable for matrix multiplication
    #tmp_sum = zero(eltype(c))
    #@time @inbounds for k = 1:size(a)[2]
    #    c[i,j] += a[i,k] * b[k, j]
    #end
    c[i,j] = dot(a[i,:], b[:, j])
    #c[i,j] = tmp_sum
end

# Creating a wrapper kernel for launching with error checks
function matmul!(a, b, c)
    #if size(a)[2] != size(b)[1]
    #    println("Matrix size mismatch!")
    #    return nothing
    #end
    #if isa(a, Array)
    kernel! = matmul_kernel!(CPU(),4)
    #else
    #    kernel! = matmul_kernel!(CUDADevice(),256)
    #end
    kernel!(a, b, c, ndrange=size(c))
    return nothing
end

n = 400
a = rand(n,n)
b = rand(n, n)
c = zeros(n, n)

@time for i=1:10 matmul!(a, b, c) end;
@time matmul!(a, b, c)
# beginning CPU tests, returns event
@time for i=1:10
    LinearAlgebra.mul!(c,a,b);
end;


@kernel function copy_kernel!(A, @Const(B))
    I = @index(Global)
    @inbounds A[I] = B[I]
end

function mycopy!(A::Array, B::Array)
    @assert size(A) == size(B)
    kernel = copy_kernel!(CPU(), 8)
    kernel(A, B, ndrange=length(A))
end

A = zeros(128, 128)
B = ones(128, 128)
@time for i = 1:10 mycopy!(A, B); end;
@time kernel(A, B, ndrange=length(A))
wait(event)
@test A == B


@time @inbounds for k = 1:size(a,2)
    @time @inbounds c[i,j] = c[i,j] + a[i,k] * b[k, j]
    #@time copy(c[i,j]);
end
@time @inbounds a[i,k]

i=1
j=1 
