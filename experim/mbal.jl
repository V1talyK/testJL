using UnicodePlots, LinearAlgebra, HypothesisTests, StatsBase
rsrc = dirname(Base.source_path())
include("$rsrc/mbal_fun.jl")

pa = 100
p00 = 100

nt = 12
qi = collect(2.5:0.5:8);  qi[1:6].= 0
qp = collect(1:12);  qp[1:6].= 4
dt = 30.5
a0 = [dt/4,1,0.2]
pf, dp = sim(pa, qi, qp, p00, a0)
pf0 = copy(pf)
pf = pf.+rand(-10:10,nt)
pfr = copy(pf)
#pfr[1:6] .= NaN
pfr[7:12] .= NaN

J = sum(abs2,pf.-pf0)
J = lf(pfr,pf0)
J_MAPE = mean(abs,(pf.-pf0)./pf0)
plt = lineplot(1:nt,pf0)
    scatterplot!(plt,pf)
    println(plt)


a1 = a0.*1.0
p, dJ_dx, dPT, d2PT, H0, x_opt = adp_fun(1000, a1, nt, pa, qi, qp, p00, pfr)
mape(p,pfr)

Δx = calc_Δx(H0,x_opt,dJ_dx,d2PT,dPT,p,pfr,nt)


Cx = inv(H0)
Cy = dJ_dx'*Cx*dJ_dx

# S2p = 1/(nt-3)*JJ[end]
# S2x = S2p./dJ_dx.^2
# Δa = tp.*sqrt.(S2x)
Δp = sqrt.(sum((dPT.*Δx).^2,dims=1)[:])
vt = 1:nt
plt = lineplot(vt,p, ylim=[10,260], xlim = [0, 36])
    lineplot!(plt, vt,p.+Δp)
    lineplot!(plt, vt,p.-Δp)
    vtf = vt[.!isnan.(pfr)]
    scatterplot!(plt,vtf, pfr[.!isnan.(pfr)])
    println(plt)

ntp = 24
vtp = length(qi).+(1:ntp)
qip = qi[end]*ones(ntp); qip[12:24] .= 10
qip = vcat(collect(qi[end]:-1:2), fill(2,5), 2:13)
qip = vcat(fill(16,6),fill(24,6),fill(32,6),fill(16,6))
qpp = qp[end]*ones(ntp);

pp, dPTp = sim(pa, qip, qpp, p[end], x_opt,dPT[:,end])
Δpp = sqrt.(sum((dPTp.*vcat(Δx,Δp[end])).^2,dims=1)[:])

lineplot!(plt,vtp,pp)
    lineplot!(plt, vtp,pp.+Δpp)
    lineplot!(plt, vtp,pp.-Δpp)
    println(plt)

tocsv(rsrc,"pres.csv",pfr, p, p.-Δp, p.+Δp)
tocsv(rsrc,"presp.csv",pp, pp.-Δpp, pp.+Δpp)
tocsv(rsrc,"qq.csv",qi, qp)
tocsv(rsrc,"qqp.csv",qip, qpp)



for (k,v) in enumerate(zip(a,Δa))
    println("$(round(v[1],digits=2))±$(round(v[2],digits=2))")
end
