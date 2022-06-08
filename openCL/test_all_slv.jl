using LinearAlgebra, OpenCL, SparseArrays, BenchmarkTools, LoopVectorization, BenchmarkTools, UnicodePlots
rsrc=dirname(dirname(Base.source_path()));
include(joinpath(rsrc,"openCL/get_Ab.jl"))
include(joinpath(rsrc,"solv/tsty.jl"))

device, ctx, queue = cl.create_compute_context()
p = cl.Program(ctx, source=low_slv_kernel) |> cl.build!
krnL = cl.Kernel(p, "slv_lowM")

p = cl.Program(ctx, source=up_slv_kernel) |> cl.build!
krnU = cl.Kernel(p, "slv_upM")

knL = make_order(U)
knU = make_order1(L)
zz = zeros(maximum(length.(knL)))

clL, rwL, nzL = Int32.(copy(U.colptr)), Int32.(copy(U.rowval)), Float32.(copy(U.nzval));
clU, rwU, nzU = Int32.(copy(L.colptr)), Int32.(copy(L.rowval)), Float32.(copy(L.nzval));

nu = 2^6
y32 = zeros(Float32,length(b),nu)
x32 = zeros(Float32,length(b),nu)
b2 = hcat([b.+i.-1 for i = 1:nu]...)

y0 = L\b2
x0 = U\y0

BLOCK_SIZE = 1024
lmem = cl.LocalMem(Float32, UInt32(BLOCK_SIZE));
zz_bf = cl.Buffer(Float32, ctx, (:rw,:use), hostbuf=Float32.(zz*0))
#row_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(list))
x_bf = cl.Buffer(Float32, ctx, (:rw, :use), hostbuf=x32)
y_bf = cl.Buffer(Float32, ctx, (:rw, :use), hostbuf=y32)

zr_bf= cl.Buffer(Float32, ctx, (:rw, :copy), hostbuf=zeros(Float32,length(x32)))
b_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(b2))
clL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=clL)
rwL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(rwL.-1))
nzL_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=nzL)

clU_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=clU)
rwU_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(rwU.-1))
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
lvl_lngL = Int32.(length.(knLn))
knUl = Int32.(vcat(knU...))
lvl_lngU = Int32.(length.(knUn))

iknL1 = Int32.(vcat(1,cumsum(length.(knLn)).+1)[1:end-1])
iknL2 = Int32.(cumsum(length.(knLn)))
knL32 = map(x->Int32.(x),knLn)

iknU1 = Int32.(vcat(1,cumsum(length.(knUn)).+1)[1:end-1])
iknU2 = Int32.(cumsum(length.(knUn)))
knU32 = map(x->Int32.(x),knUn)

knL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knL32[1])
knL1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knLl)
iknL1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=iknL1)
lvl_lngL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=lvl_lngL)


knU_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knU32[1])
knU1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knUl)
iknU1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=iknU1)
lvl_lngU_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=lvl_lngU)


#queue = cl.CmdQueue(ctx, :profile)
#@time cl.copy!(queue, zz_bf, zz_bf);
sdfL = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.([length(iknL1),L.n]))
cl.copy!(queue, y_bf, zr_bf);
gs = (BLOCK_SIZE*nu,1)
bs = (BLOCK_SIZE,1)
@time queue(krnL, gs, bs, zz_bf, knL1_bf, iknL1_bf, lvl_lngL_bf,
                            y_bf, b_bf, clL_bf, rwL_bf, nzL_bf, sdfL, lmem)
yy = cl.read(queue, y_bf)
    sum(abs,yy.-vec(y0))


sdfU = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.([length(iknU1),L.n]))
cl.copy!(queue, x_bf, zr_bf);
gs = (BLOCK_SIZE*nu,1)
bs = (BLOCK_SIZE,1)
@time queue(krnU, gs, bs, zz_bf, knU1_bf, iknU1_bf, lvl_lngU_bf,
                            x_bf, y_bf, clU_bf, rwU_bf, nzU_bf, sdfU, lmem)
xx = cl.read(queue, x_bf)
    sum(abs,xx.-vec(x0))
    mean(abs,xx.-vec(x0))


function slv_cl!(xx)
    #cl.copy!(queue, y_bf, zr_bf);
    #cl.copy!(queue, x_bf, zr_bf);
    queue(krnL, gs, bs, zz_bf, knL1_bf, iknL1_bf, lvl_lngL_bf,
                                y_bf, b_bf, clL_bf, rwL_bf, nzL_bf, sdfL, lmem)

    #queue(krn_U, gs, bs, zz_bf, knU1_bf, iknU1_bf, lvl_lngU_bf,
    #                            x_bf, y_bf, clU_bf, rwU_bf, nzU_bf, sdfU, lmem)
    #xx .= cl.read(queue, x_bf)
end

@btime slv_cl!($xx)
@btime $ACL\$b;
