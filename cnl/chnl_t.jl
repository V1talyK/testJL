function foo(vl, chnl)
  for v in vl
    v=v*2
    put!(chnl, v)
  end
end

vl = rand(100)
my_channel = Channel( (channel_arg) -> foo(vl, channel_arg) )

for i in my_channel
    i |> println
end


put!(my_channel, 2)
x = [similar(b) for i in 1:Threads.nthreads()]

bar(bb) = foo1(x,mA, bb)

mcnl1 = Channel(bar, spawn=true)
mcnl2 = Channel( (bb) -> foo1(x,mA, bb), spawn=true)

function foo1(x,mA,cnl)
  for v in cnl
    println(Threads.threadid())
    #x[Threads.threadid()] .= mA*v
    mul!(x[Threads.threadid()],mA,v)
  end
end

@time for i=1:100
  #put!(mcnl1,b)
  put!(mcnl2,b)
end

function foo2(cnl::Channel,zz::Array{Float64,1},
           x::Array{Float64,1},
           b::Array{Float64,1},
           cl::Array{Int64, 1},
           rw::Array{Int64, 1},
           nz::Array{Float64, 1})
  for v in cnl
    zz[v[2]] = calc_zz_k(v[1],x,b,cl,rw,nz)
  end
end

bar3(vk::Channel) = foo2(vk,zz,x,b,cl,rw,nz)

mcnl3 = Channel(bar3, spawn=true)
