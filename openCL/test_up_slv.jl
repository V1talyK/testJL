using LinearAlgebra, OpenCL, SparseArrays, BenchmarkTools, LoopVectorization, BenchmarkTools, UnicodePlots
rsrc=dirname(dirname(Base.source_path()));
include(joinpath(rsrc,"openCL/get_Ab.jl"))
include(joinpath(rsrc,"solv/tsty.jl"))

knU = make_order1(L)
zz = zeros(maximum(length.(knU)))
clU = copy(L.colptr); rwU = copy(L.rowval); nzU = copy(L.nzval);

x = zeros(length(b))
bp = b[ACL.p]
y0 = L\bp
x0 = U\y0

x = zeros(size(y0))
@time solv_krnU!(x,knU,y0,zz, clU,rwU,nzU)

device, ctx, queue = cl.create_compute_context()
BLOCK_SIZE = 512
p = cl.Program(ctx, source=uslv_kernel) |> cl.build!
krn_U = cl.Kernel(p, "uslvk")
#bnr = cl.info(p, :binaries);

lmem = cl.LocalMem(Float32, UInt32(maximum(length.(knU))));
zz_bf = cl.Buffer(Float32, ctx, (:rw,:use), hostbuf=Float32.(zz*0))
#row_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(list))
x32 = Float32.(x)
x_bf = cl.Buffer(Float32, ctx, (:rw, :use), hostbuf=x32)
x0_bf= cl.Buffer(Float32, ctx, (:rw, :copy), hostbuf=zeros(Float32,length(x32)))
y_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(y0))
cl1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(clU))
rw_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(rwU))
nz_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(nzU))

kn1 = Int32.(vcat(knU...))
ikn1 = Int32.(vcat(1,cumsum(length.(knU)).+1)[1:end-1])
ikn2 = Int32.(cumsum(length.(knU)))
kn32 = map(x->Int32.(x),knU)

kn_new = Vector(undef,0)
for (k,v) in enumerate(knU)
    for (k1,v1) = enumerate(Iterators.partition(v,BLOCK_SIZE))
        push!(kn_new,v1)
    end
end

kn1 = Int32.(vcat(kn_new...))
ikn1 = Int32.(vcat(1,cumsum(length.(kn_new)).+1)[1:end-1])
ikn2 = Int32.(cumsum(length.(kn_new)))
kn32 = map(x->Int32.(x),kn_new)

kn_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=kn32[1])
kn1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=kn1)
ikn1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=ikn1)
ikn2_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=ikn2)

#queue = cl.CmdQueue(ctx, :profile)
#@time cl.copy!(queue, zz_bf, zz_bf);
sdf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=[Int32(length(ikn1))])
cl.copy!(queue, x_bf, x0_bf);
bs = (BLOCK_SIZE,)
@time queue(krn_U, bs, bs, zz_bf, kn_bf, kn1_bf, ikn1_bf, ikn2_bf,
                            x_bf, y_bf, cl1_bf, rw_bf, nz_bf, sdf, lmem)
xx = cl.read(queue, x_bf)
    sum(abs,xx.-x0)
