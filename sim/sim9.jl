using SparseArrays
using  UnicodePlots, Term

function make9p(Tt, Pa, nw, bet;
                Ve = ones(nw)*100*100*10*0.2,
                lm = 1.0)
    Δt = 30.5;

    r = [1,1,2,2,2,3,3,4,4,4,5,5,5,5,6,6,6,7,7,8,8,8,9,9]
    c = [2,4,1,5,3,2,6,1,5,7,2,6,4,8,3,5,9,4,8,7,5,9,6,8]
    TT = 2 .*Tt[r].*Tt[c]./(Tt[r].+Tt[c])
    dTr = 2 .*(Tt[c]./(Tt[r].+Tt[c])).^2
    dTc = 2 .*(Tt[r]./(Tt[r].+Tt[c])).^2

    AA = sparse(r, c, TT)
    #dA = sparse(r, c, dT)

    A1 = zeros(nw)
    dA1 = zeros(nw)

    bi = [1,2,3,4,6,7,8,9]
    lam = [2,1,2,1,1,2,1,2].*lm
    bb = zeros(nw);
    bb[bi] .= Tt[bi].*Pa.*lam

    eVp = Ve*bet/Δt
    A1 .= sum(AA,dims=2)[:] + eVp;
    A1[bi] .= A1[bi] .+ Tt[bi].*lam

    #dA1 .= sum(dA, dims=1)[:];
    #dA1[bi] .= dA1[bi] .+ lam

    for v in 1:nw
        AA[v,v] = - A1[v];
        #dA[v,v] = - dA1[v];
    end
    dA = 1
    return AA, bb, eVp, dTc, dTr, r, c, lam, bi
end

function sim(qw, nt, AA, bb, P0, eVp, dTc, dTr, r, c, lam, bi)
    PM = zeros(size(qw))
    dP_dp0 = Vector(undef, nt)
    dP_dVp = Vector(undef, nt)
    dP_dT = Vector(undef, nt)
    p0 = copy(P0)
    dP_dVp0 = zeros(length(eVp),length(eVp))
    dP_dT0 = zeros(length(eVp),length(eVp))
    for t = 1:nt
        simt!(view(PM,:,t),AA,bb,view(qw,:,t), p0, eVp)
        dP_dp0[t] = zeros(length(eVp),length(eVp))
        dP_dVp[t] = zeros(length(eVp),length(eVp))
        dP_dT[t] = zeros(length(eVp),length(eVp))

        cacl_dP_dp0!(dP_dp0[t],AA,eVp)
        cacl_dP_dVp!(dP_dVp[t],AA,view(PM,:,t),p0,eVp,dP_dVp0)
        cacl_dP_dT!(dP_dT[t],AA,view(PM,:,t), p0, eVp, dP_dT0, dTc, dTr, Pa, r, c, lam, bi)

        p0 .= view(PM,:,t)
        dP_dVp0 .= dP_dVp[t]
        dP_dT0 .= dP_dT[t]
    end
    return PM, dP_dp0, dP_dVp, dP_dT
end

function simt!(Pt,AA,bt,qt, p0, eVp)
    Pt .= AA\(-bt.+qt .- p0.*eVp)
end

function cacl_dP_dp0!(dP_dp0,AA,eVp)
    bb = zeros(length(eVp))
    for (k,v) in enumerate(eVp)
        bb[k] = - eVp[k]
        dP_dp0[:,k] .= AA\bb
        bb[k] = 0.0;
    end
end

function cacl_dP_dVp!(dP_dVp,AA,pt,p0,eVp,dP_dVp0)
    bb = zeros(length(pt))
    dA_dVp = zeros(length(pt))
    for (k,v) in enumerate(pt)
        bb[k] = - p0[k]
        dA_dVp[k] = 1.0;
        dP_dVp[:,k] .= AA\(bb .+ dA_dVp.*pt .- dP_dVp0[:,k].*eVp)
        bb[k] = 0.0;
        dA_dVp[k] = 0.0;
    end
end

function cacl_dP_dT!(dP_dT,AA,pt,p0,eVp,dP_dT0, dTc, dTr, pa, r, c, lam, bi)
    db_dT = zeros(length(pt))
    dA_dT = zeros(length(pt))
    k1 = 0
    dP = pt[c] .- pt[r]
    AS = dTr.*dP# .+ dTr.*dP
    for (k,v) in enumerate(pt)
        if k in bi
            k1+=1
            db_dT[k] = (pt[k] - pa)*lam[k1]
        end
        for i in 1:length(r)
            if r[i]==k
                dA_dT[r[i]] += AS[i];
            end
            if c[i]==k
                dA_dT[r[i]] += AS[i];
            end
        end
        dP_dT[:,k] .= AA\(db_dT .- dA_dT.- dP_dT0[:,k].*eVp)
        db_dT[k] = 0.0;
        dA_dT[:] .= 0.0;
    end
end

function make_simt_f(AA,bt,eVp)
    x1 = zeros(length(bt))
    function simt_f(x,u)
        simt!(x1, AA, bt, u, x, eVp)
        return x1
    end
    return simt_f
end

panel(p; kw...) = begin
  p.margin[] = p.padding[] = 0  # make plots more compact
  Panel(string(p; color=true); style="hidden", fit=true, kw...)
end

function plot_P(PM, nw)
    plt = Vector(undef, Int64(ceil(nw/3)))
    for (k,v) in enumerate(Iterators.partition(1:nw,3))
        plt[k] = lineplot(PM[v[1],:], ylim = [floor(minimum(PM)),ceil(maximum(PM))],
                name = "скв. $(v[1])", ylabel = "P")
        for i in v[2:end]
            lineplot!(plt[k], PM[i,:], name = "скв. $i")
        end
    end
    grid(panel.(plt[1:3]); layout=(1, nothing)) |> print
end

function plot_P_aft_klm(P0, Pk, Pn, iw)
    plt = Vector(undef, 3)
    miP = min(minimum(P0[iw,:]),minimum(Pk[iw,:]),minimum(Pn[iw,:]))
    maP = max(maximum(P0[iw,:]),maximum(Pk[iw,:]),maximum(Pn[iw,:]))
    yL =  [floor(miP),ceil(maP)]

    plt[1] = lineplot(P0[iw,:], ylim = yL, name = "факт")
    lineplot!(plt[1], Pn[iw,:], name = "шум")

    plt[2] = lineplot(P0[iw,:], ylim = yL, name = "факт")
    lineplot!(plt[2], Pk[iw,:], name = "калман")

    plt[3] = lineplot(Pn[iw,:], ylim = yL, name = "шум")
    lineplot!(plt[3], Pk[iw,:], name = "калман")

    grid(panel.(plt[1:3]); layout=(1, nothing)) |> print
end

function pre_adp(T0, V0, PMf, Pa, bet;
                flag_v = true, flag_t = true, maxI = 150)
    Tt = copy(T0)
    mV = copy(V0)
    opT, opV = copy(T0),  copy(V0)
    #maxI = 150
    JJ_temp = Inf;
    JJ = zeros(maxI)
    optP = copy(PMf)
    P0 = Pa*ones(nw)
    for i=1:maxI
        AA, bb, eVp, dA, dT, r, c, lam, bi = make9p(Tt, Pa, nw, bet;
                                                Ve = 250/3*250/3*1*0.14*mV,
                                                lm = 0.0)
        PM, dP_dp0, dP_dVp, dP_dT = sim(qw, nt, AA, bb, P0, eVp, dA, dT, r, c, lam, bi)

        JJ[i] = mean(abs2, PM.-PMf)
        print(i,"  ",JJ[i])

        if JJ[i]<JJ_temp
            opT, opV = copy(Tt),  copy(mV)
            JJ_temp = JJ[i]
            optP .= PM
        end

        dJ_dT = zeros(nw)
        dJ_dV = zeros(nw)
        for t=1:nt
            dJ_dT .+= -2*dP_dT[t]'*(PMf[:,t].-PM[:,t])
            dJ_dV .+= -2*dP_dVp[t]'*(PMf[:,t].-PM[:,t])
        end
        if all(dJ_dT.==0) & all(dJ_dV.==0)
            break
        end

        if flag_t
            Tt .= Tt.*(1.0 .- 0.05.*(dJ_dT./maximum(abs, dJ_dT)))
        end
        if flag_v
            mV .= mV.*(1.0 .- 0.05.*(dJ_dV./maximum(abs, dJ_dV)))
        end

        #Tt .= Tt.*(1.0 .- 0.05.*sign.(dJ_dT))
        #mV .= mV.*(1.0 .- 0.05.*sign.(dJ_dV))
        println(" ", round.(dJ_dT, digits = 1))
    end
    lineplot(JJ)|>println
    return opT, opV, optP
end
