tst1_kernel = "__kernel void tst1 (__global uint *x1,
                                  __global uint *x2,
                                  __global uint *x3,
                                  __global uint *x4,
                                  __global uint *x5,
                                  __global uint *x6,
                                  __local float *localSums)
 {
  uint lc_id = get_local_id(0);
  uint gr_id = get_group_id(0);
  uint gl_id = get_global_id(0);

  uint wr_dim = get_work_dim(0);
  uint gl_size = get_global_size(0);
  uint gr_size = get_local_size(0);

  x1[gl_id] = gl_id;
  x2[gl_id] = lc_id;
  x3[gl_id] = gr_id;
  x4[gl_id] = gr_size;
  x5[gl_id] = gr_size;
  x6[gl_id] = gr_size;
    //for (uint j = 0; j<5; j++)
      //{
      //barrier(CLK_GLOBAL_MEM_FENCE);
      //barrier(CLK_LOCAL_MEM_FENCE);
      //mem_fence(CLK_LOCAL_MEM_FENCE);
      //mem_fence(CLK_GLOBAL_MEM_FENCE);
    //}
 }"
