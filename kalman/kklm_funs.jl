function run_klm(Pfact, Pa, qw, simt_f,dP_dp0, dP_dVp, dP_dT, dbet, dT)

    zk = Pfact .+ rand(-0.5:0.1:0.5,nw, nt).*10
        sid = mean(abs2,zk.-Pfact)

    simt_f = make_simt_f(AA,bb,eVp)
    p_kor = kalman_flt(zk, sid, Pa, qw, simt_f,dP_dp0, dP_dVp, dP_dT, dbet, dT)

    for iw = 1:3#nw
        plot_P_aft_klm(Pfact, p_kor, zk, iw)
    end

    v1 = round.(sum(abs2,Pfact.-p_kor,dims = 2)[:], digits = 1)
    v2 = round.(sum(abs2,Pfact.-zk,dims = 2)[:], digits = 1)
    println.(v1,"  ",v2, " ",v1.<v2)
end

function kalman_flt(zk, sid, Pa, qw, simt_f,dP_dp0, dP_dVp, dP_dT, dbet, dT)
    nw, nt = size(qw)
    Pk0 = zeros(nw, nw);
    Pk0[1:nw+1:nw*nw].=(2/3)^2
    Qk = [zeros(nw, nw) for _ in 1:nt]
    for t=1:nt
        #Qk[t] .= dP_dVp[t].*dbet
        #Qk[t] .= cov(dP_dVp[1][:,:],dP_dVp[1][:,:]')
        Qk[t] .= dP_dVp[t]*diagm(0=>dbet.^2 .*ones(9))*dP_dVp[t]' .+
                 dP_dT[t]*diagm(0=>dT.^2 .*ones(9))*dP_dT[t]'
    end

    Hk = diagm(0 => ones(nw))

    Rk = diagm(0 => sid*10*ones(nw))#(1/3)^2
    I1 = diagm(0 => ones(nw))#(1/3)^2

    p0 = Pa*ones(nw)
    p_sim = similar(qw)
    p_kor = similar(qw)

    dPx = dP_dp0

    for t=1:nt
        p_sim[:,t] = simt_f(p0,qw[:,t])
        Pkt  = dPx[t]*Pk0*dPx[t]' .+ Qk[t]
        yk = zk[:,t] .- p_sim[:,t]
        Sk = Hk*Pkt*Hk' + Rk
        Kk = Pkt*Hk'/Sk
        p_kor[:,t] = p_sim[:,t] + Kk*yk
        Pkk = (I1 - Kk*Hk)*Pkt
        Pk0 = Pkk
        # PP0 = PP[t]
        #p0 .= p_sim[:,t]
        p0 .= p_kor[:,t]
    end
    #plot_P(p_sim, nw)
    return p_kor
end
