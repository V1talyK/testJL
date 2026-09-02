using LinearAlgebra, SparseArrays
using Plots

function core_ppl(nx, nt, qq, pbnd)
    A = spzeros(nx, nx)
    bb = zeros(nx)
    Sp = 10.0
    Hp = 1.0
    dx = 1000/nx
    Vp = 1000/nx*Sp*Hp

    r = vcat(collect(1:(nx-1)), collect(2:nx))
    c = vcat(collect(2:nx), collect(1:(nx-1)))
    k = fill(10.0, nx)
    #k[50] = 0.05
    bet = fill(1e-3, nx)
    T = 2*k[r] .* k[c] ./ (k[r] .+ k[c])
    for (k, v) in enumerate(zip(r, c))
        A[v[1], v[2]] = T[k]/dx*Sp
    end
    for i = 1:nx
        A[i, i] = -sum(A[i, :]) - bet[i]*Vp
    end

    A[1, 1] -= k[1]/dx*Sp
    A[nx,nx] -= k[nx]/dx*Sp
    
    tmp = zeros(nx).+10
    ppl = zeros(nx, nt)
    for t=1:nt
        bb.=0
        bb[1] = -k[1]*pbnd[t,1]/dx*Sp
        bb[nx] = -k[nx]*pbnd[t,2]/dx*Sp
        bb[Int64(ceil(nx/2))] = qq[t]
        bb.-=bet .* tmp.*Vp
        tmp.=A\bb
        ppl[:, t] .= tmp
    end
    return ppl
end

nx = 100
nt = 120
pbnd = zeros(nt, 2)
qq = fill(1.0, nt)
pbnd[:,1] .= 8.0
pbnd[:,2] .= 10.0
pbnd[20:35,2] .= 9

ppl = core_ppl(nx, nt, qq, pbnd);

plt = plot(ppl[:,1])
plot!(plt, ppl[:,2])
plot!(plt, ppl[:,19])
plot!(plt, ppl[:,21])
plot!(plt, ppl[:,25])
plot!(plt, ppl[:,50])
plot!(plt, ppl[:,55])
plot!(plt, ppl[:,60])
plot!(plt, ppl[:,61])
plot!(plt, ppl[:,62])
plot!(plt, ppl[:,100], lw = 2)

plot(ppl[50,:])
plot!(ppl[5,:])

function qvazi_ppl(nx, nt, qq, pbnd)
    pql = zeros(nt)
    lam1 = 10*10*1/500
    lam2 = 10*10*1/500
    C1 = zeros(5, nt)
    C2 = zeros(5, nt)
    a = 10/1e-3

    C10 = zeros(5)
    C20 = zeros(5)
    FF = zeros(nt)
    FF = -sum(vcat(pbnd[1,:]', diff(pbnd, dims=1)), dims=2)[:]./2 + 1000/4*vcat(qq[1], diff(qq))/10

    S1 = a./omega.*(1.0 .-cos.(omega.*500 ./a))
    F1 = vcat(pbnd[1,1],diff(pbnd[:,1]))*S1

    S2 = a./omega.*(1.0 .-cos.(omega.*500 ./a))
    F2 = vcat(pbnd[1,2],diff(pbnd[:,2]))*S2

    S3 = a./omega.*(500 .- a./omega .* sin.(omega.*500 ./a))
    F3 = vcat(pbnd[1,1],diff(pbnd[:,1]))*S3

    Nk = 
    FF = - (F1+F2+F3)/Nk
    omega = (pi*collect(1:5)/(500/a + 500/a))
    omega2 = omega.^2
    for t = 1:nt
        for i=1:5
            C1[i,t] = (C10[i]+FF[t])/(1+omega2[i])
            C2[i,t] = (C20[i]+FF[t])/(1+omega2[i])
            C10[i] = C1[i, t]
            C20[i] = C2[i, t]
        end
    end
    fk = sin.(omega.*500/a).* sin.(omega.*500/a)
    for t = 1:nt
        tmp = sum(C1[:,t].*fk + C2[:,t].*fk)
        pql[t] = (pbnd[t,1]*lam1 + pbnd[t,2]*lam2 - qq[t])/(lam1+lam2) + tmp
    end
    return pql
end

pql = qvazi_ppl(nx, nt, qq, pbnd)
plot!(pql)
