
function two_step(qw, Pf)
    yy1, aa = mnk_step1(qw, Pf)
    println("1st step MASE: ", round(mean(abs2,Pf.-yy1),digits=3)," MAPE: ", round(mape(Pf,yy1), digits=3))
    yy2 = mnk_step2(qw, yy1, Pf, aa)
    println("2st step MASE: ", round(mean(abs2,Pf.-yy2),digits=3)," MAPE: ", round(mape(Pf,yy2), digits=3))
    return yy2, aa
end

function mnk_step1(xx, yy)
    nn1 = size(xx,1)
    nn2 = size(yy,1)
    xx = vcat(xx,ones(1,size(xx,2)))
    AA = xx*xx'
    aa = zeros(nn1+1,nn2)
    for i=1:nn2
        BB = xx*yy[i,:]
        aa[:,i] = AA\BB
    end
    yyr = aa'*xx
    return yyr, aa
end

function mnk_step2(xx, yy1, yy, aa)
    nn = size(yy1,1)
    xx = vcat(xx,ones(1,size(xx,2)))
    bb = zeros(nn-1+nn+1,nn)
    yyr = similar(yy)
    for i = 1:nn
        PT = copy(yy1)
        PT = PT[1:end .!=i,:]
        PT = vcat(PT,xx)
        AA = PT*PT'
        BB = PT*yy[i,:]
        bb[:,i] = AA\BB

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


function plot_P_lr(PM, P_lr, nw)
    plt = Vector(undef, Int64(ceil(nw)))
    for (k,v) in enumerate(1:nw)

        plt[v] = lineplot(PM[v,:], ylim = [floor(minimum(PM)),ceil(maximum(PM))],
                name = "fact", ylabel = "P", title = "скв. $(v)")
        #for i in v[2:end]
            lineplot!(plt[k], P_lr[v,:], name = "calc")
        #end
    end
    grid(panel.(plt); layout=(3, nothing)) |> print
end


function plot_P_lr_ex(PM, P_lr, P_ex, nw)
    plt = Vector(undef, Int64(ceil(nw)))
    for (k,v) in enumerate(1:nw)
        x_ex = size(PM,2) - size(P_ex, 2) +1 : size(PM,2)
        plt[v] = lineplot(PM[v,:], ylim = [floor(minimum(PM)),ceil(maximum(PM))],
                name = "fact", ylabel = "P", title = "скв. $(v)")
        #for i in v[2:end]
            lineplot!(plt[k], P_lr[v,:], name = "calc")
            lineplot!(plt[k], x_ex, P_ex[v,:], name = "exam")
        #end
    end
    grid(panel.(plt); layout=(3, nothing)) |> print
end
