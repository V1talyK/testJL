pm = NaN32.*zeros(Float32, nw, nt)
pm0 = fill(init_ppl, nw)
m1 = rand(Float32, nw, 2)
m2 = rand(Float32, 2, nw)
m3 = rand(Float32, nw) 
m4 = rand(Float32, nw)

loss(x) = mse(ppl.v, x)
dJ_dp = zeros(Float32, nw);
dp0_dm3 = zeros(Float32, nw, nw);
dp_dm3 = zeros(Float32, nw, nw);

dp0_dm4 = zeros(Float32, nw, nw);
dp_dm4 = zeros(Float32, nw, nw);

dJ_dm1 = zeros(Float32, nw*2)
dJ_dm2 = zeros(Float32, nw*2)
dJ_dm3 = zeros(Float32, nw)
dJ_dm4 = zeros(Float32, nw)

dM1_dm1 = zeros(Float32, nw, nw*2)
dM2_dm2 = zeros(Float32, 2, nw*2)
dM3_dm3 = zeros(Float32, nw, nw)
dM3_dm4 = zeros(Float32, nw, nw)

dM2_dM3 = zeros(Float32, 2, nw)
dM1_dM2 = zeros(Float32, nw, 2)


dp0_dm1 = zeros(Float32, nw, nw*2)
dp0_dm2 = zeros(Float32, nw, nw*2)


function run_step(m1, m2, m3, m4)
    dM2_dM3 .= m2
    dM1_dM2 .= m1

    dJ_dm1 = zeros(Float32, nw*2)
    dJ_dm2 = zeros(Float32, nw*2)
    dJ_dm3 = zeros(Float32, nw)
    dJ_dm4 = zeros(Float32, nw)

    pm.=NaN32;
    pm0 .= init_ppl
    dJ_dm1 .= 0f0;
    dJ_dm2 .= 0f0;
    dJ_dm3 .= 0f0;
    dJ_dm4 .= 0f0;
    dp0_dm1 .= 0f0;
    dp0_dm2 .= 0f0;
    dp0_dm3 .= 0f0;
    dp0_dm4 .= 0f0;

    for t = 1:16
        pm[:,t], M2, M3 = zsim(pm0, qw.v[:,t], m1, m2, m3, m4)
        pm[:,t].=clamp.(pm[:,t],-100f0, 100f0)
        dM3_dm3 .= 0f0
        dM3_dm4 .= 0f0
        dM3_dm2 = m3.*dp0_dm2
      
        dM2_dm2.=0f0
        dM2_dm2[1,1:nw] = M3
        dM2_dm2[2,nw+1:2*nw] = M3
        dM2_dm2 .+= m2*dM3_dm2
      
        dM2_dm1 = dM2_dM3*(m3.*dp0_dm1)

        dM1_dm1.=0f0
        for iw = 1:nw
            si = (iw-1)*2+1
            #println(si)
            dM1_dm1[iw,si:si+1].=M2

            dM3_dm3[iw, iw] = pm0[iw];
            dM3_dm3[:, iw] .+= m3.*dp0_dm3[:,iw];

            dM3_dm4[iw, iw] = 1f0;
            dM3_dm4[:, iw] .+= m3.*dp0_dm4[:,iw];
        end
        dM1_dm1.+=m1*dM2_dm1

        dp_dm4 = dM1_dM2*dM2_dM3*dM3_dm4
        dp_dm3 = dM1_dM2*dM2_dM3*dM3_dm3
        dp_dm2 = dM1_dM2*dM2_dm2
        dp_dm1 = dM1_dm1


        tmp = -2 .*(ppl.v[:,t] - pm[:,t])
        tmp[isnan.(tmp)].=0f0
        dJ_dp.=tmp
        
        dJ_dm1 += (dJ_dp'*dp_dm1)[:]
        dJ_dm2 += (dJ_dp'*dp_dm2)[:]
        dJ_dm3 += (dJ_dp'*dp_dm3)[:]
        dJ_dm4 += (dJ_dp'*dp_dm4)[:]
        
        pm0.=pm[:,t]
        dp0_dm1 .= dp_dm1;
        dp0_dm2 .= dp_dm2;
        dp0_dm3 .= dp_dm3;
        dp0_dm4 .= dp_dm4;
    end
    
    return loss(pm), [dJ_dm1, dJ_dm2, dJ_dm3, dJ_dm4]
end

function zsim(x0, qt, m1, m2, m3, m4)
    M3 = m3.*x0 .+ m4 .+ qt
    M2 = m2*M3
    M1 = m1*M2
    return M1, M2, M3       
end

JJ = zeros(100)
lr = 1e-4
for i=1:16
    JJ[i], dJ = run_step(m1, m2, m3, m4)
    mJ1 = quantile(abs.(reshape(dJ[1],nw, 2))[:], 0.8)
    mJ2 = quantile(abs.(reshape(dJ[2],nw, 2))[:], 0.8)
    mJ3 = quantile(abs.(dJ[3]), 0.8)
    mJ4 = quantile(abs.(dJ[4]), 0.8)

    mJ = max(mJ1, mJ2, mJ3, mJ4)
    l1 = clamp.(reshape(dJ[1],nw, 2), -mJ, mJ)#/mJ
    m1 .= m1 .- l1*lr

    mJ2 = median(abs.(reshape(dJ[2],2, nw)))
    l2 = clamp.(reshape(dJ[2],2, nw), -mJ, mJ)#/mJ
    m2 .= m2 .- l2*lr

    mJ3 = median(abs.(reshape(dJ[3], nw)))
    l3 = clamp.(reshape(dJ[3], nw), -mJ, mJ)#/mJ
    m3 .= m3 .- l3*lr

    mJ4 = median(abs.(reshape(dJ[4], nw)))
    l4 = clamp.(reshape(dJ[4], nw, 1), -mJ, mJ)#/mJ4
    m4 .= m4 .- l4*lr
    # m1 .= m1.*(1 .-0.025.*reshape(sign.(dJ[1]),nw, 2))
    # m2 .= m2.*(1 .-0.025.*reshape(sign.(dJ[2]),2, nw))
    # m3 .= m3.*(1 .-0.025.*reshape(sign.(dJ[3]),nw, 1))
    # m4 .= m4.*(1 .-0.025.*reshape(sign.(dJ[4]),nw, 1))
    println(i," ",JJ[i])
end

loss(pm)

optimizer = ADAM(0.01, 0.9, 0.999, 1e-8, nx)

# Dummy parameters and gradients
nx = nw*6
params = rand(nx)
grads = rand(nx)

# Update parameters
J = zeros(0)
p0 = vcat(m1[:], m2[:], m3[:], m4[:])
for i=1:100
    m1[:].=p0[1:nw*2]
    m2[:].=p0[nw*2+1:nw*4]
    m3[:].=p0[nw*4+1:nw*5]
    m4[:].=p0[nw*5+1:end]
    _JJ, dJ = run_step(m1, m2, m3, m4)
    grads = vcat(dJ...)
    p0 .= update!(optimizer, Float64.(p0), Float64.(grads))
    #println(p0)
    push!(J, loss(pm))
end


plot(J)