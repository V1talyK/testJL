mul_kernel = "
    #define BLOCK_SIZE $(BLOCK_SIZE)
    __kernel __attribute__((reqd_work_group_size(BLOCK_SIZE, BLOCK_SIZE, 1)))
    void matvec(__global const float *A,
                                    __global const float *x,
                                    uint ncols,
                                     __global float *y)
{
    size_t i = get_global_id(0);              // Global id, used as the row index
    size_t k = get_global_id(1);
    int local_idx0 = get_local_id(0);
    int local_idx1 = get_local_id(1);
    __global float const *a = &A[i*ncols];    // Pointer to the i'th row
    float sum = 0.f;                          // Accumulator for dot product
    for (size_t j = 0; j < ncols; j++) {
        sum += a[j] * x[j];
    }
    y[i] = local_idx0;
 }
 "

BLOCK_SIZE = 32;

a = rand(Float32, 100,20)
b = rand(Float32, 20)

device, ctx, queue = cl.create_compute_context()

a_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=a)
b_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=b)
c_buff = cl.Buffer(Float32, ctx, :w, 20)
ncol = cl.Buffer(UInt32, ctx, :r, 100)

p = cl.Program(ctx, source=mul_kernel) |> cl.build!
k = cl.Kernel(p, "matvec")

@time queue(k, size(a), (BLOCK_SIZE,BLOCK_SIZE), a_buff, b_buff, UInt32(20), c_buff)
@time c = a*b
r = cl.read(queue, c_buff)

if isapprox(norm(r .- a*b), zero(Float32))
 @info "Success!"
else
 @error "Norm should be 0.0f"
end
