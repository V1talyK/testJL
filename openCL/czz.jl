slv_kernel = "__kernel void slvk ( __global float *zz,
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
  //sizeof(ikn1)
  for (uint j = 0; j<5; j++)
      {
      if (gl_id <= ikn2[j]-ikn1[j])
         {
         uint trow = kn1[ikn1[j]-1+gl_id]-1;
         float s = 0.f;
         for (uint i = cl1[trow]-1; i<cl1[trow+1]-2; i++)
            {
            //uint zi = rw[i]-1;
            s+=x[rw[i]-1]*nz[i];
            }
         localSums[gl_id] = (b[trow]-s)/nz[cl1[trow+1]-1-1];
         //barrier(CLK_GLOBAL_MEM_FENCE);

         x[trow] = localSums[gl_id];
         }
       barrier(CLK_LOCAL_MEM_FENCE);
       barrier(CLK_GLOBAL_MEM_FENCE);
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
