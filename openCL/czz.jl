slv_kernel = "__kernel void slvk ( __global float *zz,
                                  __global const uint *kn1,
                                  __global const uint *ikn1,
                                  __global const uint *lvl_lng,
                                  __global float *x,
                                  __global const float *b,
                                  __global const uint *cl1,
                                  __global const uint *rw,
                                  __global const float *nz,
                                  __global const uint *sdf,
                         __local float *localSums)
 {
  uint local_id = get_local_id(0);
  uint group_size = get_local_size(0);
  uint gl_id = get_global_id(0);
  float s = 0.f;
  uint trow = 0;
  for (uint j = 0; j<sdf[0]; j++)
      {
      if (gl_id <= lvl_lng[j]-1)
         {
         trow = kn1[ikn1[j]-1+gl_id]-1;
         s = 0.f;
         for (uint i = cl1[trow]-1; i<cl1[trow+1]-2; i++)
            {
            s+=x[rw[i]-1]*nz[i];
            }
        x[trow] = (b[trow]-s)/nz[cl1[trow+1]-2];
        //x[trow] = localSums[gl_id];
        }
    barrier(CLK_GLOBAL_MEM_FENCE);
    }
 }"
 p = cl.Program(ctx, source=slv_kernel) |> cl.build!
 krn = cl.Kernel(p, "slvk")

slv_kernel1 = "__kernel void slvk1 ( __global float *zz,
                                   __global const uint *kn1,
                                   __global const uint *ikn1,
                                   __global const uint *lvl_lng,
                                   __global float *x,
                                   __global const float *b,
                                   __global const uint *cl1,
                                   __global const uint *rw,
                                   __global const float *nz,
                                   __global const uint *sdf,
                          __local float *localSums)
  {
   uint local_id = get_local_id(0);
   uint group_size = get_local_size(0);
   uint gl_id = get_global_id(0);
   float s = 0.f;
   uint trow = 0;
   for (uint j = 0; j<sdf[0]; j++)
       {
       if (lvl_lng[j]>1)
         {
          if (gl_id < lvl_lng[j])
            {
            trow = kn1[ikn1[j]-1+gl_id]-1;
            s = 0.f;
            for (uint i = cl1[trow]-1; i<cl1[trow+1]-2; i++)
               {
               s+=x[rw[i]]*nz[i];
               }
            x[trow] = (b[trow]-s)/nz[cl1[trow+1]-2];
            }
         }
     else
         {
         trow = kn1[ikn1[j]-1]-1;
         localSums[local_id] = 0.f;

         if (local_id < cl1[trow+1]-2 - (cl1[trow]-1))
             {
             uint i2 = cl1[trow]-1+local_id;
             localSums[local_id] = x[rw[i2]]*nz[i2];
             }

         for (uint stride = group_size/2; stride>0; stride /=2)
            {
             barrier(CLK_LOCAL_MEM_FENCE);
             if (local_id < stride)
               localSums[local_id] += localSums[local_id + stride];
            }
         if (local_id == 0)
            {
            x[trow] = (b[trow]-localSums[0])/nz[cl1[trow+1]-2];
            }
         }
     barrier(CLK_GLOBAL_MEM_FENCE);
     barrier(CLK_LOCAL_MEM_FENCE);
     }
  }"
  p = cl.Program(ctx, source=slv_kernel1) |> cl.build!
  krn = cl.Kernel(p, "slvk1")

uslv_kernel = "__kernel void uslvk ( __global float *zz,
                                   __global const uint *kn,
                                   __global const uint *kn1,
                                   __global const uint *ikn1,
                                   __global const uint *ikn2,
                                   __global float *x,
                                   __global const float *b,
                                   __global const uint *cl1,
                                   __global const uint *rw,
                                   __global const float *nz,
                                   __global const uint *sdf,
                          __local float *localSums)
  {
   uint local_id = get_local_id(0);
   uint group_size = get_local_size(0);
   uint gl_id = get_global_id(0);
   float s = 0.f;
   uint trow = 0;
   for (uint j = 0; j<sdf[0]; j++)
       {
       if (gl_id <= ikn2[j]-ikn1[j])
          {
          trow = kn1[ikn1[j]-1+gl_id]-1;
          s = 0.f;
          for (uint i = cl1[trow]; i<cl1[trow+1]-1; i++)
             {
             s+=x[rw[i]-1]*nz[i];
             }
         x[trow] = (b[trow]-s)/nz[cl1[trow]-1];
         //x[trow] = localSums[gl_id];
         }
     barrier(CLK_GLOBAL_MEM_FENCE);
     }
  }"
  p = cl.Program(ctx, source=uslv_kernel) |> cl.build!
  krn_U = cl.Kernel(p, "uslvk")

slv_kernel1 = "__kernel void slvk1 ( __global float *zz,
                                  __global const uint *kn,
                                  __global const uint *kn1,
                                  __global const uint *ikn1,
                                  __global const uint *ikn2,
                                  __global float *x,
                                  __global const float *b,
                                  __global const uint *cl1,
                                  __global const uint *rw,
                                  __global const float *nz,
                         __local float *localSums)
 {
  uint local_id = get_local_id(0);
  uint group_size = get_local_size(0);
  uint gl_id = get_global_id(0);

  int thid = get_local_id(0);
  // a_d[thid + blid * 16]

  //sizeof(ikn1)
  for (uint j = 0; j<10; j++)
      {
      barrier(CLK_GLOBAL_MEM_FENCE);
      barrier(CLK_LOCAL_MEM_FENCE);
      mem_fence(CLK_LOCAL_MEM_FENCE);
      mem_fence(CLK_GLOBAL_MEM_FENCE);
      if (gl_id <= ikn2[j]-ikn1[j])
         {
         uint trow = kn1[ikn1[j]-1+gl_id]-1;
         float s = 0.f;
         localSums[gl_id] = 0.f;
         for (uint i = cl1[trow]-1; i<cl1[trow+1]-2; i++)
            {
            s+=x[rw[i]-1]*nz[i];
            }
        if (j<10)
            {
            localSums[gl_id] = (b[trow]-s)/nz[cl1[trow+1]-1-1];
            }
        else
            {
            localSums[gl_id] = s;//(b[trow]-s)/nz[cl1[trow+1]-1-1];
            }
         if (j==1)
             {
             if(gl_id==1)
                 {
             zz[0] = cl1[trow]-1;
             zz[1] = cl1[trow+1]-2;
             zz[2] = s;
             zz[3] = x[rw[11]-1];
             zz[4] = x[rw[12]-1];
             zz[5] = rw[11];
             zz[6] = rw[12];
             }
             }
         barrier(CLK_GLOBAL_MEM_FENCE);
         barrier(CLK_LOCAL_MEM_FENCE);
         mem_fence(CLK_LOCAL_MEM_FENCE);
         mem_fence(CLK_GLOBAL_MEM_FENCE);

         if (j==1)
             {
             if(gl_id==3)
             {
                 zz[7] = x[5];
                 zz[8] = x[6];
             }
             }
         x[trow] = localSums[gl_id];
         }

    }
 }"


czz_kernel = "__kernel void czz ( __global float *zz,
                                  __global const uint *row,
                                  __global float *x,
                                  __global const float *b,
                                  __global const uint *cl1,
                                  __global const uint *rw,
                                  __global const float *nz,
                         __global float *partialSums,
                         __local float *localSums)
 {
  uint local_id = get_local_id(0);
  uint group_size = get_local_size(0);
  uint gl_id = get_global_id(0);
  float s = 0.f;
  for (uint i = cl1[row[gl_id]-1]-1; i<cl1[row[gl_id]+1-1]-2; i++)
     {
        uint zi = rw[i]-1;
        s+=x[zi]*nz[i];
     }

  x[row[gl_id]-1] = (b[row[gl_id]-1]-s)/nz[cl1[row[gl_id]+1-1]-1-1];

 }"

#//(b[row[gl_id]]-s)/nz[cl1[row[gl_id]+1]-1];
device, ctx, queue = cl.create_compute_context()

zz_buff = cl.Buffer(Float32, ctx, (:w, :copy), hostbuf=Float32.(zz))
row_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(list))
x_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(x))
b_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(b))
cl1_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(cl1))
rw_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(rw))
nz_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(nz))

gi0_buff = cl.Buffer(Int32, ctx, :w, length(a))

BLOCK_SIZE = 512
lmem = cl.LocalMem(Float32, UInt32(BLOCK_SIZE + 1));
p = cl.Program(ctx, source=czz_kernel) |> cl.build!
k = cl.Kernel(p, "czz")


@time queue(k, size(list), nothing, zz_buff, row_buff, x_buff, b_buff, cl1_buff, rw_buff, nz_buff, c_buff, lmem)
@time zzz = cl.read(queue, zz_buff)

function solv_cl_krn!(x::Array{Float64,1},
                   kn::Vector{Vector{Int64}},
                   b::Array{Float64,1},
                   zz::Array{Float64,1},
                   cl1,rw,nz)

    list = kn[1]
    for (k,row) in enumerate(list)
        sr1 = cl1[row+1]-1
        @inbounds zz[k] = b[row]/nz[sr1];
    end
    cp2!(x,zz,list)
    x_buff = cl.Buffer(Float32, ctx, (:rw, :copy), hostbuf=Float32.(x))
    for i = 2:length(kn)
        list = kn[i]

        row_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(list))

        #calc_zz!(zz,kn[i],x,b,cl1,rw,nz)
        queue(k, size(list), nothing, zz_buff, row_buff, x_buff, b_buff, cl1_buff, rw_buff, nz_buff, c_buff, lmem)
        #zz1 = cl.read(queue, zz_buff)
        #cp2!(x,zz1,list)
    end
        x .= cl.read(queue, x_buff)
end

solv_cl_krn!(x,kn,b,zz, cl1,rw,nz)
@btime solv_cl_krn!($x,$kn,$b,$zz, $cl1,$rw,$nz)
