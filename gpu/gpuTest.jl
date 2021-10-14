#/home/lik/julia/julia-0.7.0/bin/julia
Pkg.add("CUSPARSE")
Pkg.build("CUSPARSE")
Pkg.add("CUDArt")
Pkg.build("CUDArt")
using CUSPARSE

using Pkg
Pkg.build("CUDA")
Pkg.add("CUDAdrv")
Pkg.build("CUDAdrv")
Pkg.build("CuArrays")

using CUDAdrv, CUDAnative, CuArrays, CuArrays.CUSPARSE
using Test
using CuArrays, CuArrays.CUSPARSE
using CUDA
CUDA.versioninfo()

cuA =  CuArrays.CUSPARSE.CuSparseMatrixCSR(A);
cuB = CuArrays.CuArray(b)
x = CuArrays.CuArray(zeros(length(b)));

function kernel_vadd(a, b, c)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    c[i] = a[i] + b[i]
end

a = round.(rand(Float32, (300, 4)) * 100)
b = round.(rand(Float32, (300, 4)) * 100)
d_a = CuArrays.CuArray(a)
d_b = CuArrays.CuArray(b)
d_c = similar(d_a)  # output array

# run the kernel and fetch results
# syntax: @cuda [kwargs...] kernel(args...)
@time @cuda threads=12 kernel_vadd(d_a, d_b, d_c)

# CUDAdrv functionality: download data
# this synchronizes the device
c = Array(d_c)
@test a+b ≈ c

x=rand(400, 400);
y=rand(400, 1);
xs = cu(x)
ys = cu(y)
d_c = similar(xs)  # output array
ys = cu[1, 2, 3]
xs_cpu = collect(xs)

@time @cuda threads=120 kernel_vadd(xs, xs, d_c)
@time x+x;
@test x+x ≈ Array(d_c)


@time xs*ys;
@time xs\ys;

x2 = zeros(length(b))
CuArrays.CUSOLVER.csrlsvlu!(A,b,x2,1e-5,Int32(0),'O')
CuArrays.CUSOLVER.csrlsvlu!(iluA,A,cuB,x,1e-5,Int32(0),'O')

CuArrays.CUSOLVER.csrlsvchol!(cuA,cuB,x,tol,zero(Cint),'O')

CuArrays.CUSOLVER.cusolverSpCreate()

info = CuArrays.CUSPARSE.cusparseSolveAnalysisInfo_t[1][]
index = 'O'
CuArrays.CUSPARSE.ic0!('N','S',cuA,info,index);
CY = CuArrays.CUSPARSE.ilu0!('N',cuA,info,index);
iluA = CuArrays.CUSPARSE.ilu02(cuA,index);
icfA = CuArrays.CUSPARSE.ic02(cuA,index);

@time iluA*cuB
@time cuA*cuB
nj = iluA-cuA;
sum(nj)
1
