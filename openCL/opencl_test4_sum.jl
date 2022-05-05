sum_kernel =  "__kernel void sumGPU ( __global const float *input,
                         __global float *partialSums,
                         __local float *localSums)
 {
  uint local_id = get_local_id(0);
  uint group_size = get_local_size(0);

  // Copy from global to local memory
  localSums[local_id] = input[get_global_id(0)];

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


n = 2^10
a = rand(Float32, n)
b = rand(Float32, n)*0

a_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=a)
b_buff = cl.Buffer(Float32, ctx, (:w, :copy), hostbuf=b)

p = cl.Program(ctx, source=sum_kernel) |> cl.build!
k = cl.Kernel(p, "sumGPU")
lmem = cl.LocalMem(Float32, UInt32(n + 1));

@time queue(k, size(a), (UInt32(512),), a_buff, b_buff, lmem)
b = cl.read(queue, b_buff)
a = cl.read(queue, a_buff)
