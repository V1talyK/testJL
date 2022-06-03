using LinearAlgebra, OpenCL, SparseArrays, BenchmarkTools, LoopVectorization, BenchmarkTools, UnicodePlots
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
BLOCK_SIZE = 512
p = cl.Program(ctx, source=slv_kernel) |> cl.build!
krn = cl.Kernel(p, "slvk")
#bnr = cl.info(p, :binaries);

lmem = cl.LocalMem(Float32, UInt32(length(kn[1])));
zz_buff = cl.Buffer(Float32, ctx, (:rw,:use), hostbuf=Float32.(zz*0))
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

kn_new = Vector(undef,0)
for (k,v) in enumerate(kn)
    for (k1,v1) = enumerate(Iterators.partition(v,BLOCK_SIZE))
        push!(kn_new,v1)
    end
end

kn1 = Int32.(vcat(kn_new...))
ikn1 = Int32.(vcat(1,cumsum(length.(kn_new)).+1)[1:end-1])
ikn2 = Int32.(cumsum(length.(kn_new)))
kn32 = map(x->Int32.(x),kn_new)

kn_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=kn32[1])
kn1_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=kn1)
ikn1_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=ikn1)
ikn2_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=ikn2)

#queue = cl.CmdQueue(ctx, :profile)
#@time cl.copy!(queue, zz_buff, zz_buff);
sdf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=[Int32(length(ikn1))])
cl.copy!(queue, x_buff, x0_buff);
bs = (BLOCK_SIZE,)
@time queue(krn, bs, bs, zz_buff, kn_buff, kn1_buff, ikn1_buff, ikn2_buff,
                            x_buff, b_buff, cl1_buff, rw_buff, nz_buff, sdf, lmem)
xx = cl.read(queue, x_buff)
    sum(abs,xx.-x0)
    #for i=1:10 println(sum(x0[kn[i]].-xx[kn[i]])); end
    #for i=1:10 println(sum(x0[kn_new[i]].-xx[kn_new[i]])); end

#for i=1:20 println(sum(x0[kn_new[i]].-xxx[kn_new[i]])); end
e1r = [sum(abs,x0[kn_new[i]].-xx[kn_new[i]]) for i=1:length(kn_new)]
    println(lineplot(e1r))

@time cl.copy!(queue, x_buff, x0_buff);
zzz = cl.read(queue, zz_buff)


function slv_cl!(xx)
    queue(krn, (BLOCK_SIZE,), (BLOCK_SIZE,), zz_buff, kn_buff, kn1_buff, ikn1_buff, ikn2_buff,
                                x_buff, b_buff, cl1_buff, rw_buff, nz_buff, sdf,lmem)
    xx .= cl.read(queue, x_buff)
    return xx
end

@btime slv_cl!(xx)
