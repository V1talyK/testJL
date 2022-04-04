using LoopVectorization
ACL = cholesky(mA)
CL = make_CL_in_julia(ACL, 1)
@time CL = make_CL_in_julia(ACL, 1)

x = similar(b)
x.=ACL\b


CL32 = make_CL_in_julia32(ACL, 1)
b32 = Float32.(b)
x32 = similar(b32)

@btime $x.=$ACL\$b
@btime ldiv_cl!($x,$CL,$b)
@btime ldiv_cl!($x32,$CL32,$b32)


@profiler for i=1:100
    ldiv_cl!(x32,CL32,b32)
end

sum(abs,x.-x32)
