zz_buff = cl.Buffer(Float32, ctx, (:w,:copy), hostbuf=Float32.(zz))
#row_buff = cl.Buffer(Int32, ctx, (:r, :copy), hostbuf=Int32.(list))
x_buff = cl.Buffer(Float32, ctx, (:rw), hostbuf=Float32.(x))
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
bnr = cl.info(p, :binaries);

queue = cl.CmdQueue(ctx, :profile)
@time cl.copy!(queue, zz_buff, zz_buff);

@time queue(k, size(kn[1]), (nothing), zz_buff, kn_buff, kn1_buff, ikn1_buff, ikn2_buff,
                            x_buff, b_buff, cl1_buff, rw_buff, nz_buff, lmem)
@time xxx = cl.read(queue, x_buff)

function slv_cl(kn)
    queue(k, size(kn[1]), (nothing), zz_buff, kn_buff, kn1_buff, ikn1_buff, ikn2_buff,
                                x_buff, b_buff, cl1_buff, rw_buff, nz_buff, c_buff, lmem)
    xx = cl.read(queue, x_buff)
    return xx
end

@btime slv_cl($kn)