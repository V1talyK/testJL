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
    bb = zeros(nw);
    bb[bi] .= Tt[bi]*Pa

    eVp = Ve*bet
    A1 .= sum(AA,dims=2)[:] + eVp;
    A1[bi] .= A1[bi] .+ Tt[bi]
    for v in 1:nw
        AA[v,v] = -A1[v];
    end

    return AA, bb, eVp
end

function sim(qw, nt, AA, bb, P0, eVp)
    PM = zeros(size(qw))
    p0 = copy(P0)
    for t = 1:nt
        simt!(view(PM,:,t),AA,bb,view(qw,:,t), p0, eVp)
        p0 .= view(PM,:,t)
    end
    return PM
end

function simt!(Pt,AA,bt,qt, p0, eVp)
    Pt .= AA\(-bt.+qt .- p0.*eVp)
end

panel(p; kw...) = begin
  p.margin[] = p.padding[] = 0  # make plots more compact
  Panel(string(p; color=true); style="hidden", fit=true, kw...)
end
