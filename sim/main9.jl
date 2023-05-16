using  UnicodePlots, StatsBase, Term, LinearAlgebra
using ForwardDiff, Distributions
include(joinpath(Base.source_path(),"../../sim/sim9.jl"));
include(joinpath(Base.source_path(),"../../kalman/load_file.jl"));

nw, nt  = 9, 33;
qw = ql2 .+ qi2

Tt0 = rand(nw)
Tt0 = 10. *8.64*1e-3 .*ones(nw)/((1+10)/2)

Pa = 10;
P0 = Pa*ones(nw)
bet = (3.7e-4 + 7.4e-3)/2;
AA, bb, eVp, dTc, dTr, r, c, lam, bi = make9p(Tt0, Pa, nw, bet;
                                            Ve = 250/3*250/3*1*0.14*ones(nw),
                                            lm = 0.0)

#Решаем прямую задачу
PM, dP_dp0, dP_dVp, dP_dT, d2P_dVp2 = sim(qw, nt, AA, bb, P0, eVp, dTc, dTr, r, c, lam, bi)
plot_P(PM, nw)
plot_P(ppl2, nw)

#Считаем погрешности
ΔV = eVp.*0.1
ΔPM_V, ΔPM_T = calc_ΔPM(dP_dVp, dP_dT, eVp, Tt0; ΔV = ΔV, mT = 0.1, print_flag = true)

#Натягиваем сову на глобус
oT, oV, oP = pre_adp(Tt0, 0.5.*ones(nw), PM, Pa, bet;
                    flag_v = true, flag_t = true, maxI = 50);

println.(round.(oV.*eVp, digits=2),"  ",round.(eVp, digits=2))
println.(round.(oT, digits=3),"  ",Tt0)
mean(abs2,oP.-PM)

lf(PM,oP)
mape(PM,oP)
mape(Tt0,oT)
mape(ones(nw),oV)

plot_P(oP, nw)

#Считаем погрешности адаптированных параметров
AA, bb, eVp, dTc, dTr, r, c, lam, bi, d2T = make9p(oT, Pa, nw, bet;
                                            Ve = 250/3*250/3*1*0.14*ones(nw).*oV,
                                            lm = 0.0)
oP1, dP_dp0, dP_dVp, dP_dT, d2P_dVp2 = sim(qw, nt, AA, bb, P0, eVp, dTc, dTr, r, c, lam, bi, d2T)
sum(abs,oP1.-oP)

ΔVp = calc_Δprm(PM, oP, oT, oV, dP_dT, dP_dVp, d2P_dVp2)

#Считаем погрешности
ΔPM_V, ΔPM_T = calc_ΔPM(dP_dVp, dP_dT, eVp, Tt0; ΔV = ΔVp, mT = 0.1, print_flag = true)


#Прямой расчёт с оптимизированными параметрами и проверка погрешности
AA, bb, eVp, dTc, dTr, r, c, lam, bi = make9p(oT, Pa, nw, bet;
                                            Ve = 250/3*250/3*1*0.14*ones(nw).*oV,
                                            lm = 0.0)
oP1, dP_dp0, dP_dVp, dP_dT, d2P_dVp2 = sim(qw, nt, AA, bb, P0, eVp, dTc, dTr, r, c, lam, bi)
lf(PM, oP1)
mape(PM, oP)

iw = 1
    plt = lineplot(PM[iw,:].-oP[iw,:])
    lineplot!(plt,ΔPM_V[iw,:])
    lineplot!(plt,ΔPM_T[iw,:])
    println(plt)

dbet = eVp.*0.1
    ΔT = Tt0.*0.1
    run_klm(PM, Pa, qw, simt_f,dP_dp0, dP_dVp, dP_dT, dbet, ΔT)

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
