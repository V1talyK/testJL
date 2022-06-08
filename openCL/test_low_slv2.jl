using LinearAlgebra, OpenCL, SparseArrays, BenchmarkTools,
 LoopVectorization, BenchmarkTools, UnicodePlots
device, ctx, queue = cl.create_compute_context()

rsrc=dirname(dirname(Base.source_path()));
include(joinpath(rsrc,"openCL/get_Ab.jl"))
include(joinpath(rsrc,"solv/tsty.jl"))


knL = make_order(U)
zz = zeros(maximum(length.(knL)))
clL, rwL, nzL = Int32.(copy(U.colptr)), Int32.(copy(U.rowval)), Float32.(copy(U.nzval));

nu = 2^8
y32 = zeros(Float32,length(b),nu)
y2_bf = cl.Buffer(Float32, ctx, (:rw, :use), hostbuf=y32)
b2 = hcat([b.+i.-1 for i = 1:nu]...)
@time y02 = L\b2

BLOCK_SIZE = 1024
lmem = cl.LocalMem(Float32, Int32(BLOCK_SIZE+1));
zz_bf = cl.Buffer(Float32, ctx, (:rw,:use), hostbuf=Float32.(zz*0))

#p = cl.Program(ctx, source=slv_kernel) |> cl.build!
#krn = cl.Kernel(p, "slvk")
#bnr = cl.info(p, :binaries);

zr_bf= cl.Buffer(Float32, ctx, (:rw, :copy), hostbuf=zeros(Float32,length(y32)))
b2_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(b2))
clL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=clL)
rwL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(rwL.-1))
nzL_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=nzL)

knLn = Vector(undef,0)
for (k,v) in enumerate(knL)
    for (k1,v1) = enumerate(Iterators.partition(v,BLOCK_SIZE))
        push!(knLn,v1)
    end
end

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
cl.copy!(queue, y2_bf, zr_bf);
gs = (BLOCK_SIZE*nu,1)
bs = (BLOCK_SIZE,1)
@time queue(krnL, gs, bs, zz_bf, knL1_bf, iknL1_bf, lvl_lng_bf,
                            y2_bf, b2_bf, clL_bf, rwL_bf, nzL_bf, sdf, lmem)
yy2 = cl.read(queue, y2_bf)
    sum(abs,y02.-reshape(yy2,length(b),size(b2,2)))
    #for i=1:10 println(sum(x0[kn[i]].-xx[kn[i]])); end
    #for i=1:10 println(sum(x0[kn_new[i]].-xx[kn_new[i]])); end

#for i=1:20 println(sum(x0[kn_new[i]].-xxx[kn_new[i]])); end
e1r = [sum(abs,x0[kn_new[i]].-xx[kn_new[i]]) for i=1:length(kn_new)]
    println(lineplot(e1r))

@time cl.copy!(queue, x_buff, x0_buff);
zzz = cl.read(queue, zz_bf)


function slv_cl!(yy)
    queue(krnL, bs, bs, zz_bf, knL1_bf, iknL1_bf, lvl_lng_bf,
                                y2_bf, b2_bf, clL_bf, rwL_bf, nzL_bf, sdf, lmem)
    #yy .= cl.read(queue, y2_bf)
    return yy
end

@btime slv_cl!($yy2);


sd = Vector(undef,length(knL))
for (k,kni) in enumerate(knL)
    sd[k] = []
    for trow in kni
        push!(sd[k],length(clL[trow]:clL[trow+1]-1))
    end
end

maximum(maximum.(sd))


j=773


for j=1:774
    trow = knLl[iknL1[j]]
    c1 = clL[trow]:clL[trow+1]
    println(length(c1))

end

y0[i1]*nzL[i2]


local_id = 866
