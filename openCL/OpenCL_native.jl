using OpenCL, pocl_jll
OpenCL.versioninfo()
1+1
dims = (2000000,)
a = round.(rand(Float32, dims) * 100)
b = round.(rand(Float32, dims) * 100)
c = similar(a)

d_a = CLArray(a)
d_b = CLArray(b)
d_c = CLArray(c)

function vadd(a, b, c)
    i = get_global_id()
    @inbounds c[i] = a[i] + b[i]
    return
end

function vadd0(a, b, c)
    for i in eachindex(a)
    @inbounds c[i] = a[i] + b[i]
    end
    return nothing
end

len = prod(dims)
@opencl global_size=len vadd(d_a, d_b, d_c);
@btime @opencl global_size=len vadd($d_a, $d_b, $d_c);
@btime vadd0($a, $b, $c);
