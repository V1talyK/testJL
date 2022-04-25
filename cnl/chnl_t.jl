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

mcnl1 = Channel( (bb) -> foo1(x,mA, bb), spawn=true)
mcnl2 = Channel( (bb) -> foo1(x,mA, bb), spawn=true)

function foo1(x,mA,cnl)
  for v in cnl
    println(Threads.threadid())
    x[Threads.threadid()] .= mA*v
  end
end

@time for i=1:10
  put!(mcnl1,b)
  #put!(mcnl2,b)
end
