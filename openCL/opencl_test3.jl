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

dot_kernel = "__kernel void tdot (__global const float *a,
                                  __global const float *b,
                         __global float *partialSums,
                         __local float *localSums)
 {
  uint local_id = get_local_id(0);
  uint group_size = get_local_size(0);

  // Copy from global to local memory
  localSums[local_id] = a[get_global_id(0)]*b[get_global_id(0)];

  // Loop for computing localSums : divide WorkGroup into 2 parts
  for (uint stride = group_size/2; stride>0; stride /=2)
     {
      // Waiting for each 2x2 addition into given workgroup
      barrier(CLK_LOCAL_MEM_FENCE);

      // Add elements 2 by 2 between local_id and local_id + stride
      if (local_id < stride)
        localSums[local_id] += localSums[local_id + stride];
     }

  // Write result into partialSums[nWorkGroups]
  if (local_id == 0)
    partialSums[get_group_id(0)] = localSums[0];
 }"


n = 2^16
a = rand(Float32, n)
b = rand(Float32, n)
c = rand(Float32, n)*0

device, ctx, queue = cl.create_compute_context()

a_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=a)
b_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=b)
c_buff = cl.Buffer(Float32, ctx, (:w, :copy), hostbuf=c)

gi0_buff = cl.Buffer(Int32, ctx, :w, length(a))
gi1_buff = cl.Buffer(Int32, ctx, :w, length(a))
li0_buff = cl.Buffer(Int32, ctx, :w, length(a))
li1_buff = cl.Buffer(Int32, ctx, :w, length(a))

BLOCK_SIZE = 512

p = cl.Program(ctx, source=dot_kernel) |> cl.build!
k = cl.Kernel(p, "tdot")


@time queue(k, size(a), (UInt32(512),), a_buff, b_buff, c_buff, lmem)
cc = cl.read(queue, c_buff)
#li0 = cl.read(queue, li0_buff)

println.(collect(zip(gi0, li0)))
@time dot(a,b)

method = function enqueue_dot_kernel(queue::cl.CmdQueue,
                               k::cl.Kernel,
                               a::cl.Buffer{Float32},
                               b::cl.Buffer{Float32},
                               c::cl.Buffer{Float32},
                               dims)
    lmem = cl.LocalMem(Float32, UInt32(BLOCK_SIZE + 1));
    cl.set_args!(k, a, b, c, lmem)
    cl.enqueue_kernel(queue, k, (dims,), (UInt32(512),))
end

queue = cl.CmdQueue(ctx, :profile)
method(queue, k, a_buff, b_buff, c_buff, length(a))
