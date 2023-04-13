function mb(q)
    nt = length(q)
    Pa = 10
    P = zeros(length(q))
    bet = 150
    lam = 1
    dt = 30
    P0 = Pa;
    for t=1:nt
        P[t] = (-q[t] + lam*Pa+bet/dt*P0)/(bet/dt+lam)
        P0 = P[t]
    end
    return P
end

function mb2(q, omg)

    nt = length(q)
    Pa = 10
    P = zeros(length(q))
    bet = 150
    lam = 1
    dt = 30
    R2 = -complex(0,1)/(omg*bet/dt)
    R2 = 1/(omg*bet/dt)
    for t=1:nt
        P[t] = (Pa*lam - q[t])/(1/R2+lam)
    end
    return P
end

q = rand(1:4,100)
q[1:25] .= 2
q[26:50] .= 1
q[51:75] .= 4
q[76:100] .= 3

 P = mb(q)
    lineplot(P) |> println
omg = 1
Fy0 = fft(q)

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
    plt = lineplot(x, q)
    lineplot!(plt, x, yr) |> println

qr = Vector(undef, N÷n3)
Pr = Vector(undef, N÷n3)
for i in 1:N÷n3
    qr[i] = ak[i] * cos.(2π*(i-1)/tay * x) .+ bk[i] * sin.(2π*(i-1)/tay * x)
    Pr[i] = mb2(qr[i], 2π*(i-1)/tay )
end

foo(x) = map(i->sum(getindex.(x,i)),1:100)

lineplot(foo(qr)) |> println
lineplot(foo(qr[3:3])) |> println

lineplot(Pr[2]) |> println
lineplot(Pr[32].+Pr[33]) |> println
lineplot(foo(Pr)) |> println
