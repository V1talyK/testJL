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
pf = pf.+rand(-5:5,nt)

J = sum(abs2,pf.-pf0)
J_MAPE = mean(abs,(pf.-pf0)./pf0)
plt = lineplot(1:nt,pf0)
    scatterplot!(plt,pf)
    println(plt)


a1 = a0.*1.0
p, dJ_dx, dPT, d2PT, H0, x_opt = adp_fun(100, a1, nt, pa, qi, qp, p00, pf)
mape(p,pf)

Δx = calc_Δx(H0,x_opt,dJ_dx,d2PT,dPT,p,pf,nt)


Cx = inv(H0)
Cy = dJ_dx'*Cx*dJ_dx

# S2p = 1/(nt-3)*JJ[end]
# S2x = S2p./dJ_dx.^2
# Δa = tp.*sqrt.(S2x)
Δp = sqrt.(sum((dPT.*Δx).^2,dims=1)[:])
plt = lineplot(1:12,p, ylim=[50,120], xlim = [0, 36])
    lineplot!(plt, 1:12,p.+Δp)
    lineplot!(plt, 1:12,p.-Δp)
    scatterplot!(plt,pf)
    println(plt)

ntp = 24
vtp = length(qi).+(1:ntp)
qip = qi[end]*ones(ntp); qip[12:24] .= 10
qpp = qp[end]*ones(ntp);

pp, dPTp = sim(pa, qip, qpp, p[end], a_opt)
Δpp = sqrt.(sum((dPTp.*Δx).^2,dims=1)[:])

lineplot!(plt,vtp,pp)
    lineplot!(plt, vtp,pp.+Δpp)
    lineplot!(plt, vtp,pp.-Δpp)
    println(plt)

tocsv(rsrc,"pres.csv",pf, p, p.-Δp, p.+Δp)
tocsv(rsrc,"presp.csv",pp, pp.-Δpp, pp.+Δpp)



for (k,v) in enumerate(zip(a,Δa))
    println("$(round(v[1],digits=2))±$(round(v[2],digits=2))")
end
