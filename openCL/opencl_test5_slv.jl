using LinearAlgebra, OpenCL, SparseArrays, BenchmarkTools, LoopVectorization, BenchmarkTools
rsrc=dirname(dirname(Base.source_path()));
include(joinpath(rsrc,"openCL/get_Ab.jl"))
include(joinpath(rsrc,"solv/tsty.jl"))

kn = make_order(U)
zz = zeros(maximum(length.(kn)))
cl1 = copy(U.colptr)
rw = copy(U.rowval)
nz = copy(U.nzval)

x = zeros(length(b))
@time solv_krn!(x,kn,b,zz, cl1,rw,nz)

device, ctx, queue = cl.create_compute_context()

zz_buff = cl.Buffer(Float32, ctx, (:w,:copy), hostbuf=Float32.(zz))
#row_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(list))
x32 = Float32.(x)
x_buff = cl.Buffer(Float32, ctx, (:rw, :use), hostbuf=x32)
x0_buff= cl.Buffer(Float32, ctx, (:rw, :copy), hostbuf=zeros(Float32,length(x32)))
b_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(b))
cl1_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(cl1))
rw_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(rw))
nz_buff = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(nz))

kn1 = Int32.(vcat(kn...))
ikn1 = Int32.(vcat(1,cumsum(length.(kn)).+1)[1:end-1])
ikn2 = Int32.(cumsum(length.(kn)))

kn32 = map(x->Int32.(x),kn)
kn_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=kn32[1])
kn1_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=kn1)
ikn1_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=ikn1)
ikn2_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=ikn2)

BLOCK_SIZE = 512
lmem = cl.LocalMem(Float32, UInt32(length(kn[1])));
p = cl.Program(ctx, source=slv_kernel) |> cl.build!
k = cl.Kernel(p, "slvk")
#bnr = cl.info(p, :binaries);

#queue = cl.CmdQueue(ctx, :profile)
#@time cl.copy!(queue, zz_buff, zz_buff);

@time queue(k, size(kn[1]), (nothing), zz_buff, kn_buff, kn1_buff, ikn1_buff, ikn2_buff,
                            x_buff, b_buff, cl1_buff, rw_buff, nz_buff, lmem)
@time xxx = cl.read(queue, x_buff)
    for i=1:10 println(sum(x0[kn[i]].-xxx[kn[i]])); end

@time cl.copy!(queue, x_buff, x0_buff);

function slv_cl(kn)
    queue(k, size(kn[1]), (nothing), zz_buff, kn_buff, kn1_buff, ikn1_buff, ikn2_buff,
                                x_buff, b_buff, cl1_buff, rw_buff, nz_buff, c_buff, lmem)
    xx = cl.read(queue, x_buff)
    return xx
end

@btime slv_cl($kn)

function slvkjl(gl_id,zz,kn,kn1,ikn1,ikn2,x,b,cl1,rw,nz)
    for j=1:1
        if gl_id <= ikn2[j]-ikn1[j]
           trow = kn1[ikn1[j]+gl_id-1];
           s = 0.0;
           for i = cl1[trow]:cl1[trow+1]-2
              s+=x[rw[i]]*nz[i];
           end
           x[trow] = (b[trow]-s)/nz[cl1[trow+1]-1]
       end
    end
end

xx = similar(x)
for gl_id = 1:length(kn[1])
    slvkjl(gl_id,zz,kn,kn1,ikn1,ikn2,xx,b,cl1,rw,nz)
end

for i=1:1 println(sum(x0[kn[i]].-xx[kn[i]])); end
