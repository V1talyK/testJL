using LinearAlgebra, OpenCL, SparseArrays, BenchmarkTools, LoopVectorization, BenchmarkTools, UnicodePlots
rsrc=dirname(dirname(Base.source_path()));
include(joinpath(rsrc,"openCL/get_Ab.jl"))
include(joinpath(rsrc,"solv/tsty.jl"))

device, ctx, queue = cl.create_compute_context()
p = cl.Program(ctx, source=slv_kernel) |> cl.build!
krnL = cl.Kernel(p, "slvk")

p = cl.Program(ctx, source=uslv_kernel) |> cl.build!
krnU = cl.Kernel(p, "uslvk")

y0 = L\b
x0 = U\y0
knL = make_order(U)
knU = make_order1(L)
zz = zeros(maximum(length.(knL)))

clL, rwL, nzL = Int32.(copy(U.colptr)), Int32.(copy(U.rowval)), Float32.(copy(U.nzval));
clU, rwU, nzU = Int32.(copy(L.colptr)), Int32.(copy(L.rowval)), Float32.(copy(L.nzval));

y = zeros(length(b));   y32 = Float32.(y);
x = zeros(length(b));   x32 = Float32.(x);
#bnr = cl.info(p, :binaries);

BLOCK_SIZE = 512
lmem = cl.LocalMem(Float32, UInt32(BLOCK_SIZE));
zz_bf = cl.Buffer(Float32, ctx, (:rw,:use), hostbuf=Float32.(zz*0))
#row_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(list))
x_bf = cl.Buffer(Float32, ctx, (:rw, :use), hostbuf=x32)
y_bf = cl.Buffer(Float32, ctx, (:rw, :use), hostbuf=y32)

zr_bf= cl.Buffer(Float32, ctx, (:rw, :copy), hostbuf=zeros(Float32,length(x32)))
b_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(b))
clL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=clL)
rwL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=rwL)
nzL_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=nzL)

clU_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=clU)
rwU_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=rwU)
nzU_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=nzU)

knLn = Vector(undef,0)
knUn = Vector(undef,0)
for (k,v) in enumerate(knL)
    for (k1,v1) = enumerate(Iterators.partition(v,BLOCK_SIZE))
        push!(knLn,v1)
    end
end

for (k,v) in enumerate(knU)
    for (k1,v1) = enumerate(Iterators.partition(v,BLOCK_SIZE))
        push!(knUn,v1)
    end
end

knLl = Int32.(vcat(knL...))
knUl = Int32.(vcat(knU...))

iknL1 = Int32.(vcat(1,cumsum(length.(knLn)).+1)[1:end-1])
iknL2 = Int32.(cumsum(length.(knLn)))
knL32 = map(x->Int32.(x),knLn)

iknU1 = Int32.(vcat(1,cumsum(length.(knUn)).+1)[1:end-1])
iknU2 = Int32.(cumsum(length.(knUn)))
knU32 = map(x->Int32.(x),knUn)

knL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knL32[1])
knL1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knLl)
iknL1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=iknL1)
iknL2_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=iknL2)


knU_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knU32[1])
knU1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knUl)
iknU1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=iknU1)
iknU2_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=iknU2)

#queue = cl.CmdQueue(ctx, :profile)
#@time cl.copy!(queue, zz_bf, zz_bf);
sdf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=[Int32(length(iknL1))])
cl.copy!(queue, y_bf, zr_bf);
bs = (BLOCK_SIZE,)
@time queue(krn, bs, bs, zz_bf, knL_bf, knL1_bf, iknL1_bf, iknL2_bf,
                            y_bf, b_bf, clL_bf, rwL_bf, nzL_bf, sdf, lmem)
yy = cl.read(queue, y_bf)
    sum(abs,yy.-y0)


sdf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=[Int32(length(iknU1))])
cl.copy!(queue, x_bf, zr_bf);
bs = (BLOCK_SIZE,)
@time queue(krn_U, bs, bs, zz_bf, knU_bf, knU1_bf, iknU1_bf, iknU2_bf,
                            x_bf, y_bf, clU_bf, rwU_bf, nzU_bf, sdf, lmem)
xx = cl.read(queue, x_bf)
    sum(abs,xx.-x0)


function slv_cl!(xx)
    sdf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=[Int32(length(iknL1))])
    cl.copy!(queue, y_bf, zr_bf);
    cl.copy!(queue, x_bf, zr_bf);
    bs = (BLOCK_SIZE,)
    queue(krn, bs, bs, zz_bf, knL_bf, knL1_bf, iknL1_bf, iknL2_bf,
                                y_bf, b_bf, clL_bf, rwL_bf, nzL_bf, sdf, lmem)

    sdf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=[Int32(length(iknU1))])
    queue(krn_U, bs, bs, zz_bf, knU_bf, knU1_bf, iknU1_bf, iknU2_bf,
                                x_bf, y_bf, clU_bf, rwU_bf, nzU_bf, sdf, lmem)
    xx .= cl.read(queue, x_bf)
end

@btime slv_cl!($xx)
@btime $ACL\$b;
