nt = 20
upFlag = falses(nt)
upFlag[[1,7,11,14,16,17]].=true;
nl = 3;

At = zeros(Int64,nt);
ki = findall(upFlag)
for (k,v) in enumerate(ki)
    At[v:nt] .= k
end

fe = spzeros(nt,nt);
up = zeros(Int64,nt,nt);
comb0 = CartesianIndex.(zeros(Int64,nl),zeros(Int64,nl),zeros(Int64,nl));
for r=1:nt
    for c=1:nt
        if r<=c<r+nl
            fe[r,c] = c-r
            up[r,c] = At[c]
        end
    end
end
