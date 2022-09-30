using UnicodePlots
pa = 100
p00 = 100

nt = 12
qi = collect(2.5:0.5:8);  qi[1:6].= 0
qp = collect(1:12);  qp[1:6].= 4
dt = 30.5
a0 = [dt/4,1,0.2]
a1 = a0.*1.2

pf, dJ_dx, HH = calc(pa, qi, qp, p00, pf, a0)
p, dJ_dx, HH = calc(pa, qi, qp, p00, pf, a0)
pf = p.+rand(-10:10,nt)
J = sum(abs2,pf.-p)
plt = lineplot(1:12,p)
    scatterplot!(plt,pf)
    println(plt)

tp = 2.2
JJ = zeros(50)
a = copy(a1)
for j=1:50
    p, dJ_dx, HH, dPT = calc(pa, qi, qp, p00, pf, a)
    JJ[j] = sum(abs2,pf.-p)
    mn = dJ_dx./(maximum(abs, dJ_dx))
    a .= a.*(1 .- 0.02*mn)
    #println(dJ_dx)
end
println(lineplot(JJ))

Cx = inv(HH)
Cy = dJ_dx'*Cx*dJ_dx


S2p = 1/(nt-3)*JJ[end]
S2x = S2p./dJ_dx.^2
Δa = tp.*sqrt.(S2x)
Δp = sqrt.(sum((dPT.*Δa).^2,dims=1)[:])
plt = lineplot(1:12,p, ylim=[50,120], xlim = [0, 36])
    lineplot!(plt, 1:12,p.+Δp)
    lineplot!(plt, 1:12,p.-Δp)
    scatterplot!(plt,pf)
    println(plt)

ntp = 24
vtp = length(qi).+(1:ntp)
qip = qi[end]*ones(ntp); qip[12:24] .= 10
qpp = qp[end]*ones(ntp);

pp, dPTp = sim(pa, qip, qpp, p[end], a)
Δpp = sqrt.(sum((dPTp.*Δa).^2,dims=1)[:])

lineplot!(plt,vtp,pp)
    lineplot!(plt, vtp,pp.+Δpp)
    lineplot!(plt, vtp,pp.-Δpp)
    println(plt)

function calc(pa, qi, qp, p00, pf, a)
    p0 = p00
    p = zeros(12)
    dJ_dx = zeros(3)
    HH = zeros(3,3)
    d2P = zeros(3,3)
    d2A = zeros(3,3)
    d2B = zeros(3,3)
    dP = zeros(3)
    dPT = zeros(3,12)

    dB = [-a[3]/(1+a[1]*a[3])^2, 0, -a[1]/(1+a[1]*a[3])^2]

    d2B[:,1] = [3*a[3].^2/(1+a[1]*a[3])^3, 0, (a[1]*a[3] - 1)/(1+a[1]*a[3])^3]
    d2B[:,2] = [0, 0, 0]
    d2B[:,3] = [(a[1]*a[3] - 1)/(1+a[1]*a[3])^3, 0, 3*a[1].^2/(1+a[1]*a[3])^3]

    B = 1. / (1+a[3]*a[1])

    for t = 1:12
        A = p0 + a[1]*(a[2]*qi[t]-qp[t]+a[3]*pa)
        p[t] = A*B
        p0 = p[t]
        Δp = pf[t]-p0

        dA = dP .+ [a[2]*qi[t]-qp[t]+a[3]*pa, a[1]*qi[t], a[1]*pa]
        dP .= dA.*B + A.*dB
        dPT[:,t] = dP;
        dJ_dx = dJ_dx.-2*Δp.*dP

        d2A[:,1] = [d2P[1,1], d2P[2,1] + qi[t], d2P[3,1] + pa]
        d2A[:,2] = [d2P[1,2] + qi[t], d2P[2,2], d2P[3,2]]
        d2A[:,3] = [d2P[1,3] + pa, d2P[2,3], d2P[3,3]]

        d2P[:,1] = d2A[:,1].*B .+ A.*d2B[:,1] + dA[1].*dB + dA.*dB[1]
        d2P[:,2] = d2A[:,2].*B .+ A.*d2B[:,2] + dA[2].*dB + dA.*dB[2]
        d2P[:,3] = d2A[:,3].*B .+ A.*d2B[:,3] + dA[3].*dB + dA.*dB[3]

        HH[:,1] = HH[:,1] .- 2*(-dP.*dP[1] + Δp*d2P[:,1])
        HH[:,2] = HH[:,2] .- 2*(-dP.*dP[2] + Δp*d2P[:,2])
        HH[:,3] = HH[:,3] .- 2*(-dP.*dP[3] + Δp*d2P[:,3])
    end
    return p, dJ_dx, HH, dPT
end


function sim(pa, qi, qp, p00, a)
    nt = length(qi)
    p0 = p00
    p = zeros(nt)
    dP = zeros(3)
    dPT = zeros(3,nt)

    dB = [-a[3]/(1+a[1]*a[3])^2, 0, -a[1]/(1+a[1]*a[3])^2]

    B = 1. / (1+a[3]*a[1])

    for t = 1:nt
        A = p0 + a[1]*(a[2]*qi[t]-qp[t]+a[3]*pa)
        p[t] = A*B
        p0 = p[t]

        dA = dP .+ [a[2]*qi[t]-qp[t]+a[3]*pa, a[1]*qi[t], a[1]*pa]
        dP .= dA.*B + A.*dB
        dPT[:,t] = dP;
    end
    return p, dPT
end




for (k,v) in enumerate(zip(a,Δa))
    println("$(round(v[1],digits=2))±$(round(v[2],digits=2))")
end
