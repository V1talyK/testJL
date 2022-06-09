using LinearAlgebra, OpenCL, SparseArrays, BenchmarkTools, UnicodePlots
device, ctx, queue = cl.create_compute_context()

rsrc=dirname(dirname(Base.source_path()));
include(joinpath(rsrc,"openCL/get_Ab.jl"))
include(joinpath(rsrc,"openCL/libs.jl"))
include(joinpath(rsrc,"solv/tsty.jl"))

knL = make_order(U)
zz = zeros(maximum(length.(knL)))
clL, rwL, nzL = Int32.(copy(U.colptr)), Int32.(copy(U.rowval)), Float32.(copy(U.nzval));

y32 = zeros(Float32,length(b))
y_bf = cl.Buffer(Float32, ctx, (:rw, :use), hostbuf=y32)
y0 = L\b

BLOCK_SIZE = 128
lmem_min_sz = get_max_wide(knL,clL)[1]
lmem = cl.LocalMem(Float32, Int32(BLOCK_SIZE+1));
zz_bf = cl.Buffer(Float32, ctx, (:rw,:use), hostbuf=Float32.(zz*0))

#p = cl.Program(ctx, source=slv_kernel) |> cl.build!
#krn = cl.Kernel(p, "slvk")
#bnr = cl.info(p, :binaries);

zr_bf= cl.Buffer(Float32, ctx, (:rw, :copy), hostbuf=zeros(Float32,length(y32)))
b_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(b))
clL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=clL)
rwL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(rwL.-1))
nzL_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=nzL)

knLn = cut_lvl_by_BS(knL,BLOCK_SIZE);

knLl = Int32.(vcat(knLn...))
lvl_lng = Int32.(length.(knLn))
iknL1 = Int32.(vcat(1,cumsum(length.(knLn)).+1)[1:end-1])
iknL2 = Int32.(cumsum(length.(knLn)))
knL32 = map(x->Int32.(x),knLn)

knL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knL32[1])
knL1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=knLl)
iknL1_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=iknL1)
lvl_lng_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=lvl_lng)

sdf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.([length(iknL1),L.n]))
cl.copy!(queue, y_bf, zr_bf);
bs = (BLOCK_SIZE,)
@time queue(krnL, bs, bs, zz_bf, knL1_bf, iknL1_bf, lvl_lng_bf,
                            y_bf, b_bf, clL_bf, rwL_bf, nzL_bf, sdf, lmem)
yy = cl.read(queue, y_bf)
    sum(abs,yy.-y0)
    #for i=1:10 println(sum(x0[kn[i]].-xx[kn[i]])); end
    #for i=1:10 println(sum(x0[kn_new[i]].-xx[kn_new[i]])); end

#for i=1:20 println(sum(x0[kn_new[i]].-xxx[kn_new[i]])); end
e1r = [sum(abs,y0[knLn[i]].-yy[knLn[i]]) for i=1:length(knLn)]
    println(lineplot(e1r))

@time cl.copy!(queue, x_buff, x0_buff);
zzz = cl.read(queue, zz_bf)


function slv_cl!(yy)
    queue(krnL, bs, bs, zz_bf, knL1_bf, iknL1_bf, lvl_lng_bf,
                                y_bf, b_bf, clL_bf, rwL_bf, nzL_bf, sdf, lmem)
    yy = cl.read(queue, y_bf)
    return yy
end

@btime slv_cl!($yy)

get_max_wide(knLn,clL)


j=773

len = zeros(Int64,length(knLn))
for j=1:length(knLn)
    trow = knLl[iknL1[j]]
    c1 = clL[trow]:clL[trow+1]
    len[j] = length(c1)
end
maximum(len)

y0[i1]*nzL[i2]


local_id = 866
