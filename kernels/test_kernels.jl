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
    local_i = @index(Local)
    I = @index(Global)
    IGl = @index(Group, Linear)

    tb_sum = @localmem T TBSize
    #while i <= length(A)
    @inbounds tb_sum[local_i] += A[I] * B[I]
      #i += TBSize * DotBlocks
    #end

    M = @uniform @groupsize()

    @print " G" IG " L" local_i " I" I " |$M"
    # @print " L" IL
    # @print " I" I
    C[local_i] = C[local_i] + A[I]*B[I]
end

function mydot!(A::Array, B::Array, C::Array)
    @assert size(A) == size(B)
    kernel = dot_kernel!(CPU(), 4)
    kernel(A, B, C, ndrange=length(A))
    return nothing
end

a = collect(0.5:0.5:10)
b = collect(1:20)
c = zeros(4)


a = rand(100)
b = rand(100)
c = zeros(4)

a = rand(100_000)
b = rand(100_000)
c = zeros(4)

@time mydot!(a,b,c)

@time dot(a,b)
sum(c)


@kernel function dot_by_i_kernel!(A, B, C)
    IG = @index(Group)
    IL = @index(Local)
    I = @index(Global)
    IGl = @index(Group, Linear)
    #@inbounds A[I]*B[I]
    @print " G" IG
    @print " L" IL
    @print " I" I
    C[IL] = C[IL] + A[I]*B[I]
end

function mydot!(A::Array, B::Array, C::Array)
    @assert size(A) == size(B)
    kernel = dot_kernel!(CPU(), 4)
    kernel(A, B, C, ndrange=length(A))
    return nothing
end




@kernel function cssk_kernel!(@Const(X), @Const(NZ), @Const(RW), S, @Const(rng))
    IG = @index(Group)
    local_i = @index(Local)
    I = @index(Global)
    IGl = @index(Group, Linear)

    tb_sum = @localmem Float64 4
    @inbounds tb_sum[local_i] = 0.0

    v = rng[I]
    z = RW[v]
    #@print v
    tb_sum[local_i] += X[z]*NZ[v]
    @synchronize
    S[local_i] = tb_sum[local_i]
end

function cssk(cl::Array, rw::Array, nz::Array, x::Array, row::Int64)
    S = zeros(Float64,4)
    sru = cl[row]
    rng = collect(sru:cl[row+1]-2)
    kernel = cssk_kernel!(CPU(), 4)
    ev = kernel(x, nz, rw, S, rng, ndrange=length(rng))
    wait(ev)
    return sum(S)
end

nn = 10000
mm = 1000
X = 1 .+ zeros(mm)
NZ = rand(nn)
RW = rand(1:mm,nn)
S = zeros(4)
rng = collect(1:mm)

@time for i=1:10000
    cssk(X, NZ, RW, S, rng)
    #S.=0
end

cssk(X, NZ, RW, S, rng)
sum(S)

@time for i=1:10000
    css2(rng,RW,NZ,X)
end

@btime css2($rng,$RW,$NZ,$X);




function css2(rng,rw::Array{Int64,1},nz::Array{Float64,1},
              x::Array{Float64,1})
    s = zero(Float64)
    @inbounds @fastmath for v in rng
        z = getindex(rw,v)
        s += x[z]*nz[v]
    end
    return s
end

for list in kn
    for (k,v) in enumerate(list)
        zz[k] = calc_zz_k(v,x,b,cl,rw,nz)
    end
end

for (k,v) in enumerate(kn[2])
    zz[k] = calc_zz_k(v,x,b,cl,rw,nz)
end
imax = argmax(zz[1:length(kn[2])])
zz[imax]

@time calc_zz_kern(kernel,kn[2],zz,x,b,cl,rw,nz)

@kernel function calc_zz_kernel!(list::Array, zz::Array, x::Array, b::Array, cl::Array, rw::Array, nz::Array)
    grp_i = @index(Group)
    lcl_i = @index(Local)
    glb_i = @index(Global)

    s = zero(eltype(zz))
    row = Int64(list[glb_i])
    s = css1(cl,rw,nz,x,row)

    sr1 = cl[row+1]-1
    s = (b[row]-s)/nz[sr1];
    zz[glb_i] = s
    @synchronize

end

kernel = calc_zz_kernel!(CPU(), 4)

function calc_zz_kern(kernel,list::Array, zz::Array, x::Array, b::Array, cl::Array, rw::Array, nz::Array)
    ev = kernel(list,zz,x,b,cl,rw,nz, ndrange=length(list))
    wait(ev)
    return nothing
end
