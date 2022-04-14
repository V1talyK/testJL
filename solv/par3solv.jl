using SparseArrays, BenchmarkTools, LinearAlgebra, LoopVectorization, TimerOutputs

const to = TimerOutput()
reset_timer!(to)


r = vcat(collect(1:9),4:9,8:9)
c = vcat(1:9,1,1:5,5:-1:4)
L = sparse(r,c,1)
f = collect(0:8)
x0 = L\f


A = sprand(100,100,0.1)
A = A+A';
A[CartesianIndex.(1:100,1:100)].=-sum(A,dims=1)[:].-1
A = -A
ACL = cholesky(mA)
L = sparse(ACL.L)
dropzeros!(L)
b = rand(100)
x0 = L\b
U = copy(L')
kn = make_order(U)
zz = zeros(maximum(length.(kn)))

cl = copy(U.colptr)
rw = copy(U.rowval)
nz = copy(U.nzval)

x = zeros(length(b))
@time solv_krn!(x,kn,b,zz, cl,rw,nz)
@btime solv_krn!($x,$kn,$b,$zz,$cl,$rw,$nz)
@btime $x0.=$L\$b
@profiler for i=1:10 solv_krn!(x,kn,b,zz,cl,rw,nz); end;
@profile solv_krn!(x,kn,b,L,U,zz, fg,cl,rw,nz)

sum(abs,x.-x0)
reset_timer!(to1::TimerOutput)


fg(y,y1,y2,y3) = fgh!(zz,y,y1,y2,y3,x,b,cl,rw,nz)

list = kn[1]
@profiler bbr(zz,fg,list)


x2 = copy(x); x2.=0
@time forward_substit!(x2, L, b, invS)
@btime forward_substit!($x2, $L, $b, $invS)

invS = zeros(Float32,L.n)
nme(invS,L)

sum(abs,x0.-x2)
@code_warntype solv_krn!(x,kn,b,zz,cl,rw,nz)
@code_warntype fg(kn[1][1])


y = zeros(Float32,CL.L.n)
@btime nme($y,$CL32.L)

x = zeros(8)
ch = Channel{Any}(4)

ts = Threads.@spawn begin
    thr_id = Threads.threadid()
    for item in ch
        x[item] = thr_id
    end
end
