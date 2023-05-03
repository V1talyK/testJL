using  UnicodePlots, StatsBase, Term, LinearAlgebra
include(joinpath(Base.source_path(),"../../sim/sim9.jl"))

nw, nt  = 9, 100;
qw = rand(nw,nt)
qw[1,:] .= 1;  qw[1,50:60] .= 2;
qw[2,:] .= 1;
qw[3,:] .= 1;
qw[4,:] .= 1;
qw[5,:] .= -8; qw[5,25:50] .= -6;
qw[6,:] .= 1;
qw[7,:] .= 1;
qw[8,:] .= 1;   qw[8,75:80] .= 5;
qw[9,:] .= 1;

Tt = rand(nw)

Pa = 10;
P0 = Pa*ones(nw)
bet = 1e-4;

AA, bb, eVp = make9p(Tt, Pa, nw, bet)

PM, dP_dp0, dP_dVp = sim(qw, nt, AA, bb, P0, eVp)
plot_P(PM, nw)

dbet = eVp.*0.5

#grid(panel.(plt[4:6]); layout=(1, nothing)) |> print
#grid(panel.(plt[7:9]); layout=(1, nothing)) |> print
Pk0 = zeros(nw, nw);
Pk0[1:nw+1:nw*nw].=(2/3)^2
Qk = [zeros(nw, nw) for _ in 1:nt]
for t=1:nt
    #Qk[t] .= dP_dVp[t].*dbet
    #Qk[t] .= cov(dP_dVp[1][:,:],dP_dVp[1][:,:]')
    Qk[t] .= dP_dVp[t]*diagm(0=>dbet.*ones(9))*dP_dVp[t]'
end

Hk = diagm(0 => ones(nw))
zk = PM .+ rand(-0.5:0.1:0.5,nw, nt).*2
sid = mean(abs2,z.-P_0)
Rk = diagm(0 => sid*10*ones(nw))#(1/3)^2
I1 = diagm(0 => ones(nw))#(1/3)^2

p0 = Paq*ones(nw)
p_sim = similar(qw)
p_kor = similar(qw)
simt_f = make_simt_f(AA,bb,eVp)
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

plot_P(p_sim, nw)

iw = 8
    plot_P_aft_klm(PM, p_sim, zk, iw)

plt = lineplot(PM[1,:], ylim = [floor(minimum(PM[1,:])), ceil(maximum(PM[1,:]))])
    lineplot!(plt, p_sim[1,:])
    println(plt)

plt = lineplot(PM[1,:], ylim = [floor(minimum(PM[1,:])), ceil(maximum(PM[1,:]))])
    lineplot!(plt, zk[1,:])
    println(plt)

plt = lineplot(p_sim[1,:], ylim = [floor(minimum(PM[1,:])), ceil(maximum(PM[1,:]))])
    lineplot!(plt, zk[1,:])
    println(plt)
