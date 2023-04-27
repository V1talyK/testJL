using SparseArrays
using  UnicodePlots, Term

function make9p(Tt, Pa, nw, bet)
    Δt = 30.5;
    Ve = ones(nw)*100*100*10*0.2

    r = [1,1,2,2,2,3,3,4,4,4,5,5,5,5,6,6,6,7,7,8,8,8,9,9]
    c = [2,4,1,5,3,2,6,1,5,7,2,6,4,8,3,5,9,4,8,7,5,9,6,8]
    TT = 2 .*Tt[r].*Tt[c]./(Tt[r].+Tt[c])

    AA = sparse(r,c,TT)
    A1 = zeros(nw)
    bi = [1,2,3,4,6,7,8,9]
    lam = [2,1,2,1,1,2,1,2]
    bb = zeros(nw);
    bb[bi] .= Tt[bi]*Pa

    eVp = Ve*bet
    A1 .= sum(AA,dims=2)[:] + eVp;
    A1[bi] .= A1[bi] .+ Tt[bi].*lam
    for v in 1:nw
        AA[v,v] = -A1[v];
    end

    return AA, bb, eVp
end

function sim(qw, nt, AA, bb, P0, eVp)
    PM = zeros(size(qw))
    dP_dp0 = Vector(undef, nt)
    dP_dVp = Vector(undef, nt)
    p0 = copy(P0)
    for t = 1:nt
        simt!(view(PM,:,t),AA,bb,view(qw,:,t), p0, eVp)
        dP_dp0[t] = zeros(length(eVp),length(eVp))
        dP_dVp[t] = zeros(length(eVp),length(eVp))
        cacl_dP_dp0!(dP_dp0[t],AA,eVp)
        cacl_dP_dVp!(dP_dVp,AA,view(PM,:,t),p0)
        p0 .= view(PM,:,t)
    end
    return PM, dP_dp0, dP_dVp
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

function cacl_dP_dVp!(dP_dVp,AA,pt,p0)
    bb = zeros(length(pt))
    dA_dVp = zeros(length(pt))
    for (k,v) in enumerate(pt)
        bb[k] = - p0[k]
        dA_dVp[k] = 1.0;
        dP_dVp[:,k] .= AA\(bb .- dA_dVp.*pt)
        bb[k] = 0.0;
        dA_dVp[k] = 0.0;
    end
end

function make_simt_f(AA,bt,eVp)
    x1 = zeros(length(bt))
    function simt_f(x,u)
        simt!(x1,AA,bt,u, x, eVp)
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
        plt[k] = lineplot(PM[v[1],:], ylim = [floor(minimum(PM)),ceil(maximum(PM))])
        for i in v[2:end]
            lineplot!(plt[k], PM[i,:])
        end
    end
    grid(panel.(plt[1:3]); layout=(1, nothing)) |> print
end
