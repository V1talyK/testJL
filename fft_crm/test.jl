using FFTW

function mb(q, Pa)
    nt = length(q)
    #Pa = 10
    P = zeros(length(q))
    bet = 150
    lam = 1
    dt = 30
    P0 = 10;
    for t=1:nt
        P[t] = (lam*Pa[t]+bet/dt*P0 - q[t]*0)/(bet/dt+lam)
        P0 = P[t]
    end
    return P
end

function mb2(q, Pa, omg)

    nt = length(q)
    #Pa = 10
    P = zeros(length(q))
    bet = 150
    lam = 1
    dt = 30
    R2 = -complex(0,1)/(omg*bet/dt)
    R2 = 1/(omg*bet/dt)
    for t=1:nt
        P[t] = (Pa[t]*lam - q[t]*0)/(1/R2+lam)
    end
    return P
end

q = rand(1:4,100)
q[1:25] .= 2
q[26:50] .= 1
q[51:75] .= 4
q[76:100] .= 3
q = ones(100)
q[51:100] .= 2
Pa = q*10
 P = mb(q, Pa)
    lineplot(P) |> println
omg = 1
Fy0 = fft(q)
Fy0 = fft(Pa)

N = length(q)
x=1:N
n3 = 3
    Fy = Fy0[1:N÷n3]

    ak =  2/N * real.(Fy)
    bk = -2/N * imag.(Fy)  # fft sign convention
    ak[1] = ak[1]/2
    yr = zeros(N,1)
    tay = maximum(x) - minimum(x)
    for i in 1:N÷n3
        yr .+= ak[i] * cos.(2π*(i-1)/tay * x) .+ bk[i] * sin.(2π*(i-1)/tay * x)
    end
    plt = lineplot(x, Pa)
    lineplot!(plt, x, yr) |> println

qr = Vector(undef, N÷n3)
Pr = Vector(undef, N÷n3)
Par = Vector(undef, N÷n3)
for i in 1:N÷n3
    qr[i] = ak[i] * cos.(2π*(i-1)/tay * x) .+ bk[i] * sin.(2π*(i-1)/tay * x)
    Par[i] = ak[i] * cos.(2π*(i-1)/tay * x) .+ bk[i] * sin.(2π*(i-1)/tay * x)
    Pr[i] = mb2(qr[i], Par[i], 2π*(i-1)/tay)
end

foo(x) = map(i->sum(getindex.(x,i)),1:100)

lineplot(foo(qr)) |> println
lineplot(foo(qr[2:2])) |> println

for i=1:10
    lineplot(Pr[i]) |> println
end

lineplot(Pr[32].+Pr[33]) |> println
plt = lineplot(foo(Pr))
    lineplot!(plt,P)
    println(plt)
