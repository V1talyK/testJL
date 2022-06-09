using LinearAlgebra, OpenCL, SparseArrays, BenchmarkTools,
BenchmarkTools, UnicodePlots
device, ctx, queue = cl.create_compute_context()

rsrc=dirname(dirname(Base.source_path()));
include(joinpath(rsrc,"openCL/get_Ab.jl"))
include(joinpath(rsrc,"openCL/libs.jl"))
include(joinpath(rsrc,"solv/tsty.jl"))

BS = [2^i for i in 6:7]
tt = zeros(8,length(BS))

lmem_min_sz = get_max_wide(knL,clL)[1]

for pwr = 1:8
    nu = 2^pwr
    y32 = zeros(Float32,length(b),nu)
    y2_bf = cl.Buffer(Float32, ctx, (:rw, :use), hostbuf=y32)
    b2 = hcat([b.+i.-1 for i = 1:nu]...)

    @time y02 = L\b2

    zz_bf = cl.Buffer(Float32, ctx, (:rw,:use), hostbuf=Float32.(zz*0))

    #p = cl.Program(ctx, source=slv_kernel) |> cl.build!
    #krn = cl.Kernel(p, "slvk")
    #bnr = cl.info(p, :binaries);

    zr_bf= cl.Buffer(Float32, ctx, (:rw, :copy), hostbuf=zeros(Float32,length(y32)))
    b2_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=Float32.(b2))
    clL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=clL)
    rwL_bf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(rwL.-1))
    nzL_bf = cl.Buffer(Float32, ctx, (:r, :copy), hostbuf=nzL)

    for (k2,BLOCK_SIZE) in enumerate(BS)

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

        lmem = cl.LocalMem(Float32, Int32(BLOCK_SIZE+1))
        sdf = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.([length(iknL1),L.n]))
        cl.copy!(queue, y2_bf, zr_bf);
        gs = (BLOCK_SIZE*nu,1)
        bs = (BLOCK_SIZE,1)
        @time queue(krnL, gs, bs, zz_bf, knL1_bf, iknL1_bf, lvl_lng_bf,
                                    y2_bf, b2_bf, clL_bf, rwL_bf, nzL_bf, sdf, lmem)
        yy2 = cl.read(queue, y2_bf)
            sum(abs,y02.-reshape(yy2,length(b),size(b2,2)))

        function slv_cl!(yy)
            queue(krnL, gs, bs, zz_bf, knL1_bf, iknL1_bf, lvl_lng_bf,
                                        y2_bf, b2_bf, clL_bf, rwL_bf, nzL_bf, sdf, lmem)
            #yy .= cl.read(queue, y2_bf)
            return yy
        end

        tt[pwr,k2] = @belapsed slv_cl!($yy2)
    end
end


function slv_cl!(yy)
    queue(krnL, gs, bs, zz_bf, knL1_bf, iknL1_bf, lvl_lng_bf,
                                y2_bf, b2_bf, clL_bf, rwL_bf, nzL_bf, sdf, lmem)
    #yy .= cl.read(queue, y2_bf)
    return yy
end


tt = @belapsed 1
