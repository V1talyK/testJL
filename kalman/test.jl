for iw = 1:3
    Tt1 = copy(Tt)

    Tt1[iw] = 1.6
    AA, bb, eVp, dA, dT, r, c, lam, bi = make9p(Tt, Pa, nw, bet*0.1)
    PM0, dP_dp0, dP_dVp, dP_dT = sim(qw, nt, AA, bb, P0, eVp, dA, dT, r, c, lam, bi)

    AA, bb, eVp, dA, dT, r, c, lam, bi = make9p(Tt1, Pa, nw, bet*0.1)
    PM1, dP_dp0, dP_dVp, dP_dT1 = sim(qw, nt, AA, bb, P0, eVp, dA, dT, r, c, lam, bi)

    for t=1:10
        println(iw," ",t," ",all(sign.(dP_dT[t][:,iw]).==sign.(PM1[:,t] .- PM0[:,t])))
    end
end

VB = Vector(undef, 0)
VB1 = Vector(undef, 0)
iw = 2
for TT = 0.1:0.1:2
    Tt1 = copy(Tt)

    Tt1[iw] = TT
    AA, bb, eVp, dA, dT, r, c, lam, bi = make9p(Tt, Pa, nw, bet*0.1)
    PM0, dP_dp0, dP_dVp, dP_dT = sim(qw, nt, AA, bb, P0, eVp, dA, dT, r, c, lam, bi)

    AA, bb, eVp, dA, dT, r, c, lam, bi = make9p(Tt1, Pa, nw, bet*0.1)
    PM1, dP_dp0, dP_dVp, dP_dT1 = sim(qw, nt, AA, bb, P0, eVp, dA, dT, r, c, lam, bi)

    for t=1:10
        println(iw," ",TT," ",t," ",all(sign.(dP_dT[t][:,iw]).==sign.(PM1[:,t] .- PM0[:,t])))
    end
    push!(VB,dP_dT1)
    push!(VB1,PM1)
end

lineplot(0.1:0.1:2, getindex.(getindex.(VB,2),2,iw)) |> println
lineplot(0.1:0.1:2, getindex.(VB1,iw,2)) |> println
lineplot(0.1:0.1:2, getindex.(VB1,5,2)) |> println

lineplot(0.1:0.1:2, getindex.(getindex.(VB,3),2,iw)) |> println
