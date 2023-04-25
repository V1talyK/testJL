using SparseArrays


nw, nt  = 9, 100;
Δt = 30.5;
Tt = rand(nw)
Ve = ones(nw)*100*100*10*0.2
Pa = 10;
P0 = Pa*ones(nw)
bet = 1e-5;
qw = rand(nw,nt)
qw[5,:] = -5*qw[5,:]
r = [1,1,2,2,2,3,3,4,4,4,5,5,5,5,6,6,6,7,7,8,8,8,9,9]
c = [2,4,1,5,3,2,6,1,5,7,2,6,4,8,3,5,9,4,8,7,5,9,6,8]
TT = 2 .*Tt[r].*Tt[c]./(Tt[r].+Tt[c])
AA = sparse(r,c,TT)
bi = [1,2,3,4,6,7,8,9]
ba = zeros(nw); ba[bi] .= Tt[bi]*Pa
for v in bi
    AA[v,v] = -Tt[v] - sum(AA[v,:])
end

AA\(-ba.-qw[:,1])
