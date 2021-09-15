bpp = copy(bp);
forward_substitution!(y1, LL, bp)
x1 = backward_substitution!(Lt, y1)

for i=1:16610
    LL[i,i] = 0
end

sum(abs,y.-y1)
sum(abs,x.-x1[p])


LL*y1

sum(LL[4,:].*y1)
sum(LL[4,:].*y)


sum(LL[3,:].*y1) - bp[3]
sum(LL[3,:].*y) - bp[3]


for n = 100:100:1000
    MM = sprandn(n,n,0.25)
    MM[1:MM.n+1:MM.n^2] .= sqrt(n)
    LL = sparse(LowerTriangular(MM))
    UU = sparse(UpperTriangular(MM))

    b = rand(n)
    y1 = similar(b)
    x1 = similar(b)
    frw_sb!(y1, LL, b)
    brw_sb!(x1, UU, y1)

    y0 = LL\b
    x0 = MM\b

    println(sum(abs,x1.-x0),"     ",sum(abs,y1.-y0))
end

@time frw_sb!(x1,S, b)
@time frw_sbt!(x1,S, b)
@time x0 = S\b


@btime frw_sb!($x1,$S, $b)
@btime frw_sbt!($x1,$S, $b)
@btime $x0 .= $S\$b

@profiler for i = 1: 100
    frw_sb!(x1, S, b)
end

@btime x0 = ACL\b;

x1 = zeros(length(x0))
x2 = zeros(length(x0))
LL = sparse(ACL.L)
UU = sparse(LL')
bp = b[ACL.p]
@btime begin
    #@inbounds x2p = view(x2,ACL.p)
    forward_substit!(x1,LL, bp)
    backward_substit!(x2,UU, x1)
end


bp32 = Float32.(bp)
x132 = Float32.(x1)
x232 = Float32.(x2)
@btime begin
    forward_substit!(x132,LL32, bp32)
    backward_substit!(x232,UU32, x132)
end

sum(abs,x0.-x1)

@btime $x0 .= $ACL\$b
@btime ldiv!($x0,$ALU,$b)


@btime ldiv_cl!($x1,$CL,$b)
@profiler for i = 1:10 ldiv_cl!(x1,CL,b); end
@profiler for i = 1:100 forward_substit!(x132,LL32, bp32); end
tt = zeros(Int64,183)
xx = zeros(Float32,183)
@btime forward_substit!(x132,LL32, bp32, tt, xx);


@btime forward_substit!(x1,LL, bp)
@btime backward_substit!(x2,UU, x1)


@profiler CL = make_CL_in_julia(ACL);
@time cholesky(mA)
