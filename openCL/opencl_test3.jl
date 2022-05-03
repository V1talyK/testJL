function enqueue_naive_kernel{T}(queue::cl.CmdQueue,
                                 k::cl.Kernel,
                                 dst::cl.Buffer{T},
                                 src::cl.Buffer{T},
                                 dims)
    h, w = dims
    @assert w % BLOCK_SIZE == 0
    @assert h % BLOCK_SIZE == 0
    cl.set_args!(k, dst, src, uint32(h), uint32(h))
    cl.enqueue_kernel(queue, k, (h, w), (BLOCK_SIZE, BLOCK_SIZE))
end

naive_transpose = "
    #define BLOCK_SIZE $(BLOCK_SIZE)
    __kernel
    __attribute__((reqd_work_group_size(BLOCK_SIZE, BLOCK_SIZE, 1)))
    void transpose(__global *a_t,
                        __global *a,
                        unsigned a_width,
                        unsigned a_height)
   {
   int read_idx = get_global_id(0) + get_global_id(1) * a_width;
   int write_idx = get_global_id(1) + get_global_id(0) * a_height;
   a_t[write_idx] = a[read_idx];
   }";

dot_kernel = "
   #define BLOCK_SIZE $(BLOCK_SIZE)
    __kernel
    __attribute__((reqd_work_group_size(BLOCK_SIZE, 1, 1)))
     void tdot(__global const float *a,
                __global const float *b,
                __global int *gi0,
                __global int *li0,
                __local float *a_local)
{
    size_t gl_id = get_global_id(0);              // Global id, used as the row index
    size_t lc_id = get_local_id(0);
    size_t gr_id = get_group_id();

    gi0[gl_id] = gl_id;
    li0[gl_id] = lc_id;

    int base_idx_a   = get_group_id(0) * BLOCK_SIZE + get_group_id(1) * (BLOCK_SIZE * 1);

    int glob_idx_a   = base_idx_a + get_local_id(0) + 1 * get_local_id(1);

    a_local[lc_id] += a[gl_id]+b[gl_id];

    barrier(CLK_LOCAL_MEM_FENCE);

    gi0[0] = a_local[lc_id];
 }
 "


a = rand(Float32, 32)
b = rand(Float32, 32)

device, ctx, queue = cl.create_compute_context()

a_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=a)
b_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=b)
gi0_buff = cl.Buffer(Int32, ctx, :w, length(a))
gi1_buff = cl.Buffer(Int32, ctx, :w, length(a))
li0_buff = cl.Buffer(Int32, ctx, :w, length(a))
li1_buff = cl.Buffer(Int32, ctx, :w, length(a))

p = cl.Program(ctx, source=dot_kernel) |> cl.build!
k = cl.Kernel(p, "tdot")
lmem = cl.LocalMem(Float32, UInt32(BLOCK_SIZE + 1));

@time queue(k, size(a), (UInt32(16),), a_buff, b_buff, gi0_buff, li0_buff, lmem)
gi0 = cl.read(queue, gi0_buff)
li0 = cl.read(queue, li0_buff)

println.(collect(zip(gi0, li0)))
