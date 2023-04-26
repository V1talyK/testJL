using UnicodePlots, StatsBase, Term
include(joinpath(Base.source_path(),"../../mb/mb.jl"))
panel(p; kw...) = begin
  p.margin[] = p.padding[] = 0  # make plots more compact
  Panel(string(p; color=true); style="hidden", fit=true, kw...)
end

nt = 100
q = ones(nt); q[Int32(nt/2):end] .= 0.5
Paq = 10;
db = 100
dl = 0.1
P_0, dPx, dPu, dP_db, dP_dl, mbt_f = mb(q,Paq; lam = 0.5)
Qk = sqrt.((dP_db.*db).^2 .+ (dP_dl.*dl).^2)

Pk = zeros(length(q))
P = zeros(length(q))
PP = zeros(length(q))
Pk0 = (2/3)^2
Q = (2/3)^2
z = P_0 .+ rand(-0.5:0.1:0.5,nt)
sid = mean(abs2,z.-P_0)
z[50:60] .= 10
Rk = sid*10#(1/3)^2

PP0 = Paq
for t=1:nt
    P[t] = mbt_f(PP0,q[t])
    Pk[t] = dPx[t]*Pk0*dPx[t] + Q
    yk = z[t] - P[t]
    Sk = 1*Pk[t]*1 + Rk
    Kk = Pk[t]*1/Sk
    PP[t] = P[t] + Kk*yk
    Pkk = (1 - Kk*1)*Pk[t]
    Pk0 = Pkk
    PP0 = PP[t]
end

plt1 = lineplot(P_0,ylim = [7,10], title = "Фильтрованный")
       lineplot!(plt1,1:nt, P);

plt2 = lineplot(P_0,ylim = [7,10], title = "Исходный")
       lineplot!(plt2,1:nt, z);

grid(panel.([plt2, plt1]); layout=(1, nothing)) |> print

sum(abs2,z.-P_0)
sum(abs2,P_0.-P)

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
    P_opt0 = P_opt[t]
end
plt = lineplot(Ki)
    println(plt)

plt = lineplot(P_opt; name ="klm")
    lineplot!(plt,P; name = "dflt")
    scatterplot!(plt,P1; name = "rnd")
    println(plt)
