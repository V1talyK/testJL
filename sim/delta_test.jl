SPM = Vector(undef, 100)
for ii = 1:100
    if ii==1
        Tt = Tt0.- ΔT
    elseif ii == 2
        Tt = Tt0.+ ΔT
    else
        Tt = Tt0.*(1 .+ ΔT.*rand(-1:1,9))
    end

    AA, bb, eVp, dA, dT, r, c, lam, bi = make9p(Tt, Pa, nw, bet;
                                            Ve = 250/3*250/3*1*0.14*ones(nw),
                                            lm = 0.0)

    PM, dP_dp0, dP_dVp, dP_dT = sim(qw, nt, AA, bb, P0, eVp, dA, dT, r, c, lam, bi)
    SPM[ii] = copy(PM)
end

DPM = Vector(undef, 100)
for ii = 1:100
    DPM[ii] = abs.(SPM[ii] .- PM)
end

df = Vector(undef, nw)
for iw = 1:nw
    df[iw] = maximum(hcat(getindex.(DPM,iw,:)...),dims=2)[:]
    #df[iw] = mean(hcat(getindex.(DPM,iw,:)...),dims=2)[:]
    plt = lineplot(df[iw])
        lineplot!(plt, ΔPM[iw,:])
        plt |> println
    mean(abs, (ΔPM[iw,:] .- df[iw])./ΔPM[iw,:]) |> println
    cor(ΔPM[iw,:],df[iw]) |> println
end



function calc_Δx(H0,x0,dJ_dx,d2PT,dPT,p,pf,nt)
    x = H0\(H0*x0.-dJ_dx)
    d2P_dx2 = d2PT'
    dH_dp = LinearAlgebra.diagm(0 => -2*sum(d2P_dx2,dims=1)[:])
    dJx_dp = -2*sum(dPT,dims=2)[:]

    N = count(.!isnan.(pf))
    Jf = sum(abs2,filter(!isnan,p.-pf))/(N-9)

    dx_dp = H0\(dH_dp*(x0-x)-dJx_dp)
    Sx = sqrt.(dx_dp.^2 .*Jf)
    Sxs = Sx./sqrt(nt)

    #tp = ifelse(N==6,2.57,2.2)
    tp = quantile(TDist(N-1),0.975)
    Δx = tp*Sxs

    for (k,v) in enumerate(zip(x0,Δx))
        println("$(round(v[1],digits=3))±$(round(v[2],digits=3))")
    end
    return Δx
end
