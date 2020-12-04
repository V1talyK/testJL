using SparseArrays, UnicodePlots
nt = 20
upFlag = falses(nt)
upFlag[[1,7,11,14,16,17]].=true;
nl = 3;

function make_calc_prog(upFlag,nl)
    nt = length(upFlag)
    At = zeros(Int64,nt);
    ki = findall(upFlag)
    for (k,v) in enumerate(ki)
        At[v:nt] .= k
    end

    fe = spzeros(Int64,nt,nt);
    up = zeros(Int64,nt,nt);
    cf = spzeros(Int64,nt,nt);
    gt = spzeros(Int64,nt,nt);
    cmb = CartesianIndex.(zeros(Int64,nl),zeros(Int64,nl),zeros(Int64,nl),zeros(Int64,nl));
    rc = Vector(undef,0)
    Ro = Vector(undef,nt); for i=1:nt Ro[i] = []; end;
    for r=1:nt
        for c=r:nt
            if r<=c<r+nl
                push!(rc,[r,c])
                fe[r,c] = c-r
                up[r,c] = At[c]
                upfe = c>1 ? up[r,c-1] : 0
                upse = c>2 ? up[r,c-2] : 0

                cf[r,c] = CartesianIndex(fe[r,c],up[r,c],upfe,upse) != cmb[c-r+1]
                cmb[c-r+1] = CartesianIndex(fe[r,c],up[r,c],upfe,upse)
                if cf[r,c]!=1
                    gt[r,c] = c-r+1
                end
                push!(Ro[c],[fe[r,c],up[r,c],cf[r,c],gt[r,c]])
            end
        end
    end

    fev = [fe[irc[1],irc[2]] for irc in rc]
    upv = [up[irc[1],irc[2]] for irc in rc]

    ci = CartesianIndex.(fev,upv)
    unique(ci)
    return Ro
end
