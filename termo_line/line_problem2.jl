function ppl_lin(n, nt, bet0, q, p0, pL)
    Lx = 1000
    dx = Lx/n
    S = 10
    dV = dx*S
    A = spzeros(n, n)
    bb = zeros(n)
   

    r = vcat(collect(1:(n-1)), collect(2:n))
    c = vcat(collect(2:n), collect(1:(n-1)))
    k = fill(100.0, n)
    #k[Int64(ceil(n/2))] = 5.5
    bet = fill(bet0, n)
    T = 2*k[r] .* k[c] ./ (k[r] .+ k[c])
    for (k, v) in enumerate(zip(r, c))
        A[v[1], v[2]] = T[k]/dx*S
    end
    for i = 1:n
        A[i, i] = -sum(A[i, :]) - bet[i]*dV
    end

    A[1, 1] -= k[1]/dx*S
    A[n, n] -= k[n]/dx*S


    #p0[50:70].=15
    tmp = zeros(n)
    tmp.=10
    ppl = zeros(n, nt)
    for t=1:nt
        bb.=0
        bb[1] = -k[1]*p0[t]/dx*S
        bb[n] = -k[n]*pL[t]/dx*S
        bb[Int64(ceil(n/2))] += q[t]
        bb.-=bet .* tmp*dV
        tmp.=A\bb
        ppl[:, t] .= tmp
    end
    return ppl
end

nt = 120
q = fill(1.0, nt)
#q[20:80] .= 0
q[80:90] .= 2

p0 = fill(10, nt)
pL = fill(10, nt)
pL[50:60].=15

n = 100
bet0 = 1e-2
ppl = ppl_lin(n, nt, bet0, q, p0, pL);
plt = plot(ppl[Int64(ceil(n/2)),:])
plot!(plt, ppl[Int64(ceil(n/2)),:])


plt = plot(ppl[:,1])
plot!(plt, ppl[:,5])
plot!(plt, ppl[:,49])
plot!(plt, ppl[:,50])
plot!(plt, ppl[:,55])
plot!(plt, ppl[:,60])
plot!(plt, ppl[:,61])
plot!(plt, ppl[:,62])
plot!(plt, ppl[:,100], lw = 2)

ppl1 = zeros(nt)
lam1 = k[1]*S/500
lam2 = k[1]*S/500
ppl0 = 10
vp = S*100*bet0
for t=1:nt
    ppl1[t] = (ppl0 - (q[t] - lam1*p0[t] - lam2*pL[t])/vp)/(1+lam1/vp+lam2/vp)
    ppl0 = ppl1[t]
end

plt = plot(ppl1)
plot!(plt, ppl1)
ic = Int64(ceil(n/2))
plot!(plt, ppl[ic,:])
plot!(plt, mean(ppl[ic-1:ic+1,:], dims=1)[:])
plot!(plt, mean(ppl, dims=1)[:])

C5 = zeros(nt, 5)
F5 = zeros(nt, 5)
C50 = zeros(5)
L = 500
l5 = pi.*(2.0.*collect(1:5).-1)./(2L)
dmu = vcat(p0[1]-init_ppl, diff(p0))
dq = vcat(q[1], diff(q))
for j=1:5
    F5[:,j] = 2.0./(L.*l5[j]).*(dmu.*(-1).^j + L./k[1].*dq)
end
a = k[1]/bet0
for t=1:nt
    for j=1:5
        C5[t,j] = (C50[j] + F5[t,j])/(1+a*l5[j]^2)
        C50[j] = C5[t,j]
    end
    
    
end

plot(ppl1.+sum(C5, dims=2)[:])
plot!(ppl1)
plot!(ppl[ic,:])