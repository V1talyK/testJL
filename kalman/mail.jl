using UnicodePlots
eVp = 1e3
Δt = 30.5
nt = 100
P = zeros(nt)
P0 = 10
q = ones(nt); q[Int32(nt/2):end] .= 0.5
Paq = 10;
lam = 1

for t=1:nt
    P[t] = (P0 - (q[t] - lam*Paq)*Δt/eVp)/(1+lam*Δt/eVp)
    P0 = P[t]
end

plt = lineplot(P,ylim = [8,12])
    println(plt)
P1 = P .+ rand(-0.5:0.1:0.5,nt)

lineplot!(plt,1:nt, P1)
    println(plt)



for t=1:nt
    T1 = 1/(1+lam*Δt/eVp)
    T2 = -T1*Δt/eVp
    T3 = -T2*lam
    P[t] = T1*P0 + T2*q[t] + T3*Paq
    P0 = P[t]
end

Ki = zeros(nt)
sig1 = 1
sig2 = 1
Eek = sig2
P_opt = zeros(nt)
P_opt0 = 10
qt = 0
for t = 1:nt
    Eek1 = sig2*(Eek+sig1)/(Eek + sig1 + sig2)
    Ki[t] = Eek1/sig2
    P_opt[t] = Ki[t]*P1[t] + (1-Ki[t])*(P_opt0 + qt)
    Eek = Eek1
    qt = q[t]
end
plt = lineplot(Ki)
    println(plt)

plt = lineplot(P_opt; name ="klm")
    lineplot!(plt,P; name = "dflt")
    lineplot!(plt,P1; name = "rmd")
    println(plt)
