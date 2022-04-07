using KernelAbstractions, Test, CUDA, LinearAlgebra, BenchmarkTools

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

n = 100
a = rand(n,n)
b = rand(n, n)
c = zeros(n, n)

@time for i=1:10 matmul!(a, b, c) end;
@time matmul!(a,b, c)
# beginning CPU tests, returns event
@btime LinearAlgebra.mul!($c,$a,$b);



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


@kernel function dot_kernel!(A, B, C)
    IG = @index(Group)
    IL = @index(Local)
    I = @index(Global)
    IGl = @index(Group, Linear)
    #@inbounds A[I]*B[I]
    #@print " G" IG
    #@print " L" IL
    #@print " I" I
    C[IL] = C[IL] + A[I]*B[I]
end

function mydot!(A::Array, B::Array, C::Array)
    @assert size(A) == size(B)
    kernel = dot_kernel!(CPU(), 4)
    event = kernel(A, B, C, ndrange=length(A))
    wait(event)
    return sum(C)
end

a = collect(0.5:0.5:8)
b = collect(1:16)
c = zeros(4)


a = rand(100)
b = rand(100)
c = zeros(4)

a = rand(10_000)
b = rand(10_000)
c = zeros(4)

@time mydot!(a,b,c)

dot(a,b)
sum(c)
