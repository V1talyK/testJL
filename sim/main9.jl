using  UnicodePlots, StatsBase, Term, LinearAlgebra
include(joinpath(Base.source_path(),"../../sim/sim9.jl"))

nw, nt  = 9, 33;
qw = ql2 .+ qi2

Tt = rand(nw)
Tt = 10. *8.64*1e-3 .*ones(nw)/((1+10)/2)

Pa = 10;
P0 = Pa*ones(nw)
bet = (3.7e-4 + 7.4e-3)/2;
AA, bb, eVp, dA, dT, r, c, lam, bi = make9p(Tt, Pa, nw, bet;
                                            Ve = 250/3*250/3*1*0.14*ones(nw),
                                            lm = 0.0)

PM, dP_dp0, dP_dVp, dP_dT = sim(qw, nt, AA, bb, P0, eVp, dA, dT, r, c, lam, bi)
plot_P(PM, nw)
plot_P(ppl2, nw)

dbet = eVp.*0.1
    dT = Tt.*0.1
    run_klm(PM, Pa, qw, simt_f,dP_dp0, dP_dVp, dP_dT, dbet, dT)

dbet = eVp.*0.9
    dT = Tt.*0.8
    run_klm(ppl2, Pa, qw, simt_f,dP_dp0, dP_dVp, dP_dT, dbet, dT)

plt = lineplot(PM[1,:], ylim = [floor(minimum(PM[1,:])), ceil(maximum(PM[1,:]))])
    lineplot!(plt, p_sim[1,:])
    println(plt)

plt = lineplot(PM[1,:], ylim = [floor(minimum(PM[1,:])), ceil(maximum(PM[1,:]))])
    lineplot!(plt, zk[1,:])
    println(plt)

plt = lineplot(p_sim[1,:], ylim = [floor(minimum(PM[1,:])), ceil(maximum(PM[1,:]))])
    lineplot!(plt, zk[1,:])
    println(plt)
