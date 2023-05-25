function test_1()
    #Тестируем, задаём отклонение проводимости в 10% и смотрим какое посчитается
    AA, bb, eVp, dTc, dTr, r, c, lam, bi, d2T = make9p(Tt0, Pa, nw, bet;
                                                Ve = 250/3*250/3*1*0.14*ones(nw),
                                                lm = 0.0)
    #Решаем прямую задачу
    PM, dP_dp0, dP_dVp, dP_dT, d2P_dVp2 = sim(qw, nt, AA, bb, P0, eVp, dTc, dTr, r, c, lam, bi, d2T)

    mn = 0.9
    oT = mn.*Tt0
    #oT = copy(Tt0); oT[1] = 0.9*oT[1];
    oV = ones(nw)
    AA, bb, eVp, dTc, dTr, r, c, lam, bi, d2T = make9p(oT, Pa, nw, bet;
                                            Ve = 250/3*250/3*1*0.14*ones(nw).*oV,
                                            lm = 0.0)
                                            oP1, dP_dp0, dP_dVp, dP_dT, d2P_dVp2, d2P_dT2 = sim(qw, nt, AA, bb, P0, eVp, dTc, dTr, r, c, lam, bi, d2T)
    println("------------")
    println(round(1 .- mn, digits=2),"  ",round(sum(abs,oP1.-PM), digits=1))
    ΔVp, ΔT = calc_Δprm(PM, oP1, oT, oV, dP_dT, dP_dVp, d2P_dT2, d2P_dVp2)
end

function test_2()
    #Тестируем, задаём отклонение параметра и смотрим какое посчитается
    AA, bb, eVp, dTc, dTr, r, c, lam, bi, d2T = make9p(Tt0, Pa, nw, bet;
                                                Ve = 250/3*250/3*1*0.14*ones(nw),
                                                lm = 0.0)
    #Решаем прямую задачу
    PM, dP_dp0, dP_dVp, dP_dT, d2P_dVp2 = sim(qw, nt, AA, bb, P0, eVp, dTc, dTr, r, c, lam, bi, d2T)

    for mn = 0.85:0.05:1.15
        oT = mn.*Tt0
        #oT = copy(Tt0); oT[1] = 0.9*oT[1];
        oV = ones(nw)
        AA, bb, eVp, dTc, dTr, r, c, lam, bi, d2T = make9p(oT, Pa, nw, bet;
                                                Ve = 250/3*250/3*1*0.14*ones(nw).*oV,
                                                lm = 0.0)
                                                oP1, dP_dp0, dP_dVp, dP_dT, d2P_dVp2, d2P_dT2 = sim(qw, nt, AA, bb, P0, eVp, dTc, dTr, r, c, lam, bi, d2T)
        println("------------")
        println(round(1 .- mn, digits=2),"  ",round(sum(abs,oP1.-PM), digits=1))
        ΔVp, ΔT = calc_Δprm(PM, oP1, oT, oV, dP_dT, dP_dVp, d2P_dT2, d2P_dVp2)
    end

end
