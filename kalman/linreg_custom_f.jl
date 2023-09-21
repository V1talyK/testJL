function two_step(qw, Pf)
    yy1, aa, Δaa, wrn_flag, rnk = mnk_step1(qw, Pf)
    println("1st step MASE: ", round(mean(abs2,Pf.-yy1),digits=3)," MAPE: ", round(mape(Pf,yy1), digits=3))
    yy2 = mnk_step2(qw, yy1, Pf, aa)
    println("2st step MASE: ", round(mean(abs2,Pf.-yy2),digits=3)," MAPE: ", round(mape(Pf,yy2), digits=3))
    return yy2, aa, Δaa, wrn_flag, rnk
end

function mnk_step1(xx, yy; critical_value = 0.95)
    nn1 = size(xx,1)
    nn2 = size(xx,2)
    nny = size(yy,1)
    xx = vcat(xx,ones(1,size(xx,2)))
    AA = xx*xx'
    rnk = rank(AA)
    println("ранг1: ",rnk)
    wrn_flag = rnk!=nn1+1
    #ev = eigvals(AA)
    #println(ev," ",extrema(ev))
    aa = zeros(nn1+1,nny)
    for i=1:nny
        BB = xx*yy[i,:]
        aa[:,i] = AA\BB
    end
    yyr = aa'*xx

    sig = sum(abs2,yyr[:] .- yy[:])/(nn2 - (nn1+1))
    #println(inv(AA))
    #println(sig)
    #println(cov(xx'))
    if !wrn_flag
        CovM3x = inv(AA)*sig
        #println.(eachrow(CovM3x))
        tp = quantile(TDist(nn2 - (nn1+1)),(1+critical_value)/2)
        #println(diag(CovM3x))
        Δaa = sqrt.(diag(CovM3x))*tp
    else
        Δaa = NaN*zeros(nn1+1)
    end

    return yyr, aa, Δaa, wrn_flag, rnk
end

function mnk_step2(xx, yy1, yy, aa)
    nx = size(xx,1)
    nn = size(yy1,1)
    xx = vcat(xx,ones(1,size(xx,2)))
    bb = zeros(nx-1+nn+1,nn)
    yyr = similar(yy)
    for i = 1:nn
        PT = copy(yy1)
        PT = PT[1:end .!=i,:]
        PT = vcat(PT,xx)
        AA = PT*PT'
        BB = PT*yy[i,:]
        bb[:,i] = AA\BB
        println("ранг2: ",rank(AA))
        yyr[i,:] = bb[:,i]'*PT #+ aa[:,i]'*xx
    end
    return yyr
end

function mnk_step3(xx, yy2, Pf, aa)
    V = cov(Pf.-yy2, dims = 2)
    println(V)
    # nn = size(yy1,1)
    # xx = vcat(xx,ones(1,size(xx,2)))
    # bb = zeros(nn-1,nn)
    # yyr = similar(yy)
    # for i = 1:nn
    #     PT = copy(yy1)
    #     PT = PT[1:end .!=i,:]
    #     #PT[i,:] .= 1.0
    #     AA = PT*PT'
    #     BB = PT*(yy[i,:].-(aa[:,i]'*xx)')
    #     bb[:,i] = AA\BB
    #
    #     yyr[i,:] = bb[:,i]'*PT + aa[:,i]'*xx
    # end
    return yyr
end

function calc_rank_of_data(xx, fact_i)
    n = size(xx,1)
    MR = zeros(n+1,n+1)
    MR[1:n,1:n].= view(xx,:,fact_i)*view(xx,:,fact_i)'
    MR[n+1,1:n].=sum(view(xx,:,fact_i),dims=2)
    MR[1:n,n+1].=MR[n+1,1:n]
    MR[n+1,n+1] = 1
    return rank(MR)
end
