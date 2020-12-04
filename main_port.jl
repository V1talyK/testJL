function getPort(x)
    i_port = findfirst(x.=="--port")
    val = x[i_port+1]
end

println(ARGS)
val = getPort(ARGS)
println(val)




#fun, IN, OUT = getInputCmd(ARGS);
