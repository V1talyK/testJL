using  UnicodePlots, StatsBase, Term, LinearAlgebra
using ForwardDiff, Distributions
include(joinpath(Base.source_path(),"../../sim/sim9.jl"));
include(joinpath(Base.source_path(),"../../kalman/load_file.jl"));

nw, nt  = 9, 33;
qw = ql2 .+ qi2
qw .= mean(qw, dims = 2)
qw[qw.>0] .= mean(qw[qw.>0])
qw[qw.<0] .= mean(qw[qw.<0])

Tt0 = rand(nw)
Tt0 = 10. *8.64*1e-3 .*ones(nw)/((1+10)/2)

Pa = 10;
P0 = Pa*ones(nw)
bet = (3.7e-4 + 7.4e-3)/2;
AA, bb, eVp, dTc, dTr, r, c, lam, bi, d2T = make9p(Tt0, Pa, nw, bet;
                                            Ve = 250/3*250/3*1*0.14*ones(nw),
                                            lm = 0.0)

#Решаем прямую задачу
PM, dP_dp0, dP_dVp, dP_dT, d2P_dVp2 = sim(qw, nt, AA, bb, P0, eVp, dTc, dTr, r, c, lam, bi, d2T)
plot_P(PM, nw)
plot_P(ppl2, nw)

#Считаем погрешности
ΔV = eVp.*0.1
ΔT = Tt0.*0.1
ΔPM_V, ΔPM_T = calc_ΔPM(dP_dVp, dP_dT, eVp, Tt0; ΔV = ΔV, ΔT = ΔT, print_flag = true)

#Натягиваем сову на глобус
oT, oV, oP = pre_adp(0.9.*Tt0, ones(nw), PM, Pa, bet;
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
oT = 0.9.*Tt0
oV = ones(nw)
AA, bb, eVp, dTc, dTr, r, c, lam, bi, d2T = make9p(oT, Pa, nw, bet;
                                            Ve = 250/3*250/3*1*0.14*ones(nw).*oV,
                                            lm = 0.0)
oP1, dP_dp0, dP_dVp, dP_dT, d2P_dVp2, d2P_dT2 = sim(qw, nt, AA, bb, P0, eVp, dTc, dTr, r, c, lam, bi, d2T)
sum(abs,oP1.-oP)

ΔVp, ΔT = calc_Δprm(PM, oP, oT, oV, dP_dT, dP_dVp, d2P_dT2, d2P_dVp2)

#Считаем погрешности
ΔPM_V, ΔPM_T = calc_ΔPM(dP_dVp, dP_dT, eVp, Tt0; ΔV = ΔVp, ΔT = ΔT, print_flag = true)


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

du_dp = zeros(nw,length(PM))
inx = collect(Iterators.product(1:nw,1:nt))[:]
LI = LinearIndices((1:nw,1:nt))
for iw = 1:nw
    for iw2 = 1:nw
        for t = 1:nt
            ip = LI[iw2, t]
            z1 = 0.0
            z2 = 0.0
            for t1 = 1:nt
                for ii = 1:nw
                    for jj = 1:nw
                        z1 += dP_dT[t1][ii,iw]*dP_dT[t1][ii,jj]
                        z2 += (PM[ii,t1].-oP1[ii,t1])*d2P_dT2[t1][jj][ii,iw]
                    end
                end
            end
            du_dp[iw, ip] = sum(dP_dT[t][iw2,:])/(-z1 - z2)
        end
    end
end

N = count(.!isnan.(PM))
Jf = sum(abs2,filter(!isnan,PM.-oP1))/(N-nw)

Sx = sqrt.(sum(du_dp.^2, dims=2)[:])
Sxs = sqrt(Jf).*Sx

#tp = ifelse(N==6,2.57,2.2)
tp = quantile(TDist(N-9),0.975)
Δx = tp*Sxs
