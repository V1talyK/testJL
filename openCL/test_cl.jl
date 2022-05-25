using LinearAlgebra, OpenCL, SparseArrays, BenchmarkTools, LoopVectorization, BenchmarkTools

n = 100
device, ctx, queue = cl.create_compute_context()
x1, x2, x3, x4, x5, x6 = [zeros(Int32,n) for i in 1:6]

x1_buff = cl.Buffer(Int32, ctx, (:w,:copy), hostbuf=x1)
x2_buff = cl.Buffer(Int32, ctx, (:w,:copy), hostbuf=x2)
x3_buff = cl.Buffer(Int32, ctx, (:w,:copy), hostbuf=x3)
x4_buff = cl.Buffer(Int32, ctx, (:w,:copy), hostbuf=x4)
x5_buff = cl.Buffer(Int32, ctx, (:w,:copy), hostbuf=x5)
x6_buff = cl.Buffer(Int32, ctx, (:w,:copy), hostbuf=x6)

BLOCK_SIZE = 16
lmem = cl.LocalMem(Float32, UInt32(4));
p = cl.Program(ctx, source=tst1_kernel) |> cl.build!
k = cl.Kernel(p, "tst1")


queue(k, (n,), (UInt32(10),), x1_buff, x2_buff, x3_buff, x4_buff, x5_buff, x6_buff, lmem)
x_1 = cl.read(queue, x1_buff)
x_2 = cl.read(queue, x2_buff)
x_3 = cl.read(queue, x3_buff)
x_4 = cl.read(queue, x4_buff)
x_6 = cl.read(queue, x6_buff)
