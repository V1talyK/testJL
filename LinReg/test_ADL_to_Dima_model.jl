ql1
ppl1

qq = ql1.+qi1
cc = 1 .- qo1./qq
cc = 1 .- cumsum(qo1,dims=2)./cumsum(qq, dims = 2)
tgap = 600
iw = 2

xx = vcat(qq[[1,2,3,4,5,6,7,8,9],:],cc[[2],:],ones(1,982))
xx = vcat(qq[[1,2,3,4,5,6,7,8,9],:],qw1[[2,4,5,6,8],:],qo1[[2,4,5,6,8],:],ones(1,982))

xxa = xx[:,1:tgap]
xxe = xx[:,tgap:end]

AA = xxa*xxa'
bb = xxa*ppl1[iw,1:tgap]

aa = AA\bb

pplc = xx'*aa
mean(abs,(ppl1[iw,:].-pplc)./ppl1[iw,:])

plt = lineplot(ppl1[iw,:], ylim = [0, 20])
    lineplot!(plt, pplc)
    println(plt)

Δppl = diff(ppl1[iw,:])

p_ind = [1,2,3,4,5,6,7,8,9]
xxa1 = hcat(ppl1[p_ind,1:tgap-1]',xxa[:,2:tgap]')
xxe1 = hcat(ppl1[p_ind,tgap:end-1]',xx[:,tgap+1:end]')

AA1 = xxa1'*xxa1
bb1 = xxa1'*Δppl[1:tgap-1]

aa1 = AA1\bb1

Δpplc_a = xxa1*aa1
Δpplc_e = xxe1*aa1

pplc[2:tgap] = ppl1[iw,1] .+ cumsum(Δpplc_a)
pplc[tgap+1:end] = pplc[tgap] .+ cumsum(Δpplc_e)

using MultivariateStats

M = fit(PCA, xxa; maxoutdim = 20)
Y_fact = predict(M, xxa)
Y_miss = predict(M, xxe)


Y_fact*Δppl[1:tgap-1]

sd = (Y_fact*Y_fact')\(Y_fact*Δppl[1:tgap])

(sd'*Y_miss)'
pplc[tgap:end] = pplc[tgap] .+ cumsum((sd'*Y_miss)')
