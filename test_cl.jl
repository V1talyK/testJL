using OpenCL
import OpenCL.cl.CLArray
import CLBLAS
const clblas = CLBLAS


clblas.setup()
device, ctx, queue = clblas.get_next_compute_context()
alpha = 1.0
beta = 0.0
queue =1

hA = rand(Float32,5, 10)
hB = rand(10, 5)
A = CLArray(hA)
B = CLArray(queue, hB)
C = cl.zeros(queue, 5, 5)

clblas.gemm!('N', 'N', alpha, A, B, beta, C)
hC = cl.to_host(C)
if isapprox(hC, hA * hB)
    info("Success!")
else
    error("Results diverged")
end

using CLArrays

for dev in CLArrays.devices()
    CLArrays.init(dev)
    x = zeros(CLArray{Float32}, 5, 5) # create a CLArray on device `dev`
end
