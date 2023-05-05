Pa = 10;
P0 = Pa*ones(nw)

#Тест_1
Tt = rand(0.01:0.01:1, nw)
bet = (3.7e-4 + 7.4e-3)/2;

AA, bb, eVp, dA, dT, r, c, lam, bi = make9p(Tt, Pa, nw, bet;
                                            Ve = 250/3*250/3*1*0.14*ones(nw),
                                            lm = 0.0)

PM, dP_dp0, dP_dVp, dP_dT = sim(qw, nt, AA, bb, P0, eVp, dA, dT, r, c, lam, bi)
plot_P(PM, nw)

oT, oV, oP = pre_adp(0.5*ones(nw), ones(9), PM; flag_v = false);
println.(round.(oT, digits=1),"  ",Tt)
mean(abs2,oP.-PM)
plot_P(oP, nw)


#Тест_2
Tt = rand(0.01:0.01:1, nw)
mV = ones(nw)
bet = (3.7e-4 + 7.4e-3)/2;

AA, bb, eVp, dA, dT, r, c, lam, bi = make9p(Tt, Pa, nw, bet;
                                            Ve = 250/3*250/3*1*0.14*mV,
                                            lm = 0.0)

PM, dP_dp0, dP_dVp, dP_dT = sim(qw, nt, AA, bb, P0, eVp, dA, dT, r, c, lam, bi)
plot_P(PM, nw)

oT, oV, oP = pre_adp(Tt, 0.5*ones(9), PM; flag_v = true, flag_t = false);
println.(round.(oV, digits=1),"  ",mV)
mean(abs2,oP.-PM)
plot_P(oP, nw)
