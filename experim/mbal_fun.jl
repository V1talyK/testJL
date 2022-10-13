
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
    d2PT = zeros(3,12)

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
        Δp = ifelse(isnan(Δp),0.,Δp)

        dA = dP .+ [a[2]*qi[t]-qp[t]+a[3]*pa, a[1]*qi[t], a[1]*pa]
        dP .= dA.*B + A.*dB
        dPT[:,t] = dP;
        dJ_dx .= dJ_dx.-2*Δp.*dP

        d2A[:,1] = [d2P[1,1],         d2P[2,1] + qi[t], d2P[3,1] + pa]
        d2A[:,2] = [d2P[1,2] + qi[t], d2P[2,2],         d2P[3,2]]
        d2A[:,3] = [d2P[1,3] + pa,    d2P[2,3],         d2P[3,3]]

        d2P[:,1] = d2A[:,1].*B .+ A.*d2B[:,1] + dA[1].*dB + dA.*dB[1]
        d2P[:,2] = d2A[:,2].*B .+ A.*d2B[:,2] + dA[2].*dB + dA.*dB[2]
        d2P[:,3] = d2A[:,3].*B .+ A.*d2B[:,3] + dA[3].*dB + dA.*dB[3]
        # println(d2A[:,2])
        # println(d2A[2,:])
        # println(B)
        # println("-----------")
        d2PT[:,t] .= d2P[[1,5,9]];

        HH[:,1] = HH[:,1] .- 2*(-dP.*dP[1] + Δp*d2P[:,1])
        HH[:,2] = HH[:,2] .- 2*(-dP.*dP[2] + Δp*d2P[:,2])
        HH[:,3] = HH[:,3] .- 2*(-dP.*dP[3] + Δp*d2P[:,3])
    end
    return p, dJ_dx, HH, dPT, d2PT
end


function sim(pa, qi, qp, p00, a,dP0)
    nt = length(qi)
    p0 = p00
    p = Vector(undef,nt)
    dP = zeros(4)
    dP[1:3].=dP0
    dP[4] = 1.
    dA = zeros(4)
    dPT = zeros(4,nt)
    dPT = Vector(undef,nt)

    B = 1. / (1+a[3]*a[1])
    dB = [-a[3]/(1+a[1]*a[3])^2, 0., -a[1]/(1+a[1]*a[3])^2, 0.]

    for t = 1:nt
        A = p0 + a[1]*(a[2]*qi[t]-qp[t]+a[3]*pa)
        p[t] = A*B
        p0 = p[t]

        dA = dP .+ [a[2]*qi[t]-qp[t]+a[3]*pa, a[1]*qi[t], a[1]*pa, 0.]
        dP = dA.*B + A.*dB
        dPT[t] = dP;
    end
    return p, dPT
end

function adp_fun(maxI, a0, nt, pa, qi, qp, p00, pf0)
    nx = length(a0)
    JJ = zeros(maxI)
    a = copy(a0)
    a_opt = copy(a0)
    dJ_dx = zeros(length(a))
    dPT = zeros(nx,nt)
    d2PT = zeros(nx,nt)
    p = zeros(nt)
    HH = zeros(nx, nx)
    max_dJ_dx = 0.
    for j=1:maxI
        p, dJ_dx, HH, dPT, d2PT = calc(pa, qi, qp, p00, pf0, a)
        JJ[j] = sum(abs2,filter(!isnan,pf0.-p))
        mn = dJ_dx./max(maximum(abs, dJ_dx),max_dJ_dx)
        max_dJ_dx = max(maximum(abs, dJ_dx),max_dJ_dx)
        a_opt .= a;
        a .= a.*(1 .- 0.10*mn)
        #println(dJ_dx," ",a_opt)
    end
    println(lineplot(JJ))
    return  p, dJ_dx, dPT, d2PT, HH, a_opt
end

function tocsv(rsrc,nm,v...)
    outfile = string(rsrc,"/$nm")
    open(outfile, "w") do f
    nline = length(v[1])
      for i in 1:nline
          for v1 in v
              print(f, round(v1[i],digits = 3), "\t")
          end
          println(f,"")
      end
    end # the fi
end

function calc_Δx(H0,x0,dJ_dx,d2PT,dPT,p,pf,nt)
    x = H0\(H0*x0.-dJ_dx)
    d2P_dx2 = d2PT'
    dH_dp = LinearAlgebra.diagm(0 => -2*sum(d2P_dx2,dims=1)[:])
    dJx_dp = -2*sum(dPT,dims=2)[:]

    N = count(.!isnan.(pf))
    Jf = sum(abs2,filter(!isnan,p.-pf))/(N-3)

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

mape(xf,xc) = mean(abs,filter(!isnan,(xf.-xc)./xf))
lf(xf,xc) = sum(abs2,filter(!isnan,xf.-xc))
