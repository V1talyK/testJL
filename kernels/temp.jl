for list in kn
    for (k,v) in enumerate(list)
        zz[k] = calc_zz_k(v,x,b,cl,rw,nz)
    end
end

for (k,v) in enumerate(kn[12])
    zz[k] = calc_zz_k(v,x,b,cl,rw,nz)
end
imax = argmax(zz[1:length(kn[10])])
zz[imax]

zz.=0

const TBSize = 1024::Int
const DotBlocks = 256::Int
T = Float64

@kernel function kdot(@Const(a), @Const(b), sz, partial)
    local_i = @index(Local)
    group_i = @index(Group)
    tb_sum = @localmem T TBSize
    @inbounds tb_sum[local_i] = 0.0

    # do dot first
    i = @index(Global)
    while i <= sz
      @inbounds tb_sum[local_i] += a[i] * b[i]
      i += TBSize * DotBlocks
    end

    # then tree reduction
    offset = @private Int64 (1,)
    @inbounds begin
      offset[1] = @groupsize()[1] ÷ 2
      while offset[1] > 0
        @synchronize
        if (local_i - 1) < offset[1]
          tb_sum[local_i] += tb_sum[local_i+offset[1]]
        end
        offset[1] ÷= 2
      end
    end

    if (local_i == 1)
      @inbounds partial[group_i] = tb_sum[local_i]
    end
  end

function mydot!(a::Array, b::Array, sz, partial)
    @assert size(a) == size(b)
    kernel = kdot(CPU(), 4)
    ev = kernel(a, b, sz ,partial, ndrange=length(a))
    wait(ev)
    return nothing
end


a = [1,2,3,4]
b = [1,2,3,1]
partial = zeros(1)
sz = 4
mydot!(a::Array, b::Array, sz, partial)
# driver
#wait(kdot(CPU(), 4)(a, b, size, partial_sum, ndrange = TBSize * DotBlocks))
