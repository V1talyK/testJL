Pkg.add("CUSPARSE")
Pkg.add("CUDArt")
using CUSPARSE

Pkg.build("CUDAnative")
Pkg.build("CUDAdrv")
Pkg.add("CuArrays")

using CUDAdrv, CUDAnative, CuArrays
using Test

function kernel_vadd(a, b, c)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    c[i] = a[i] + b[i]
end

a = round.(rand(Float32, (300, 4)) * 100)
b = round.(rand(Float32, (300, 4)) * 100)
d_a = CuArray(a)
d_b = CuArray(b)
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
@time xs.\ys;
