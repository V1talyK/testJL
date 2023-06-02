PM
qw

a = 2
b = -5
x = rand(0:0.1:10,100)
y = a.*x .+ b
y1 = y.+rand(length(y))

plt = lineplot(x,y)
    scatterplot!(plt,x,y1)
    println(plt)

[100 sum(x); sum(x) sum(x.^2)]\[sum(y1), sum(x.*y1)]

a = [1, 3, 4]
x = rand(0:0.1:10,100,3)
y = x*a
y1 = y.+rand(length(y))

AA = x'*x
BB = x'*y1

AA\BB

yy1, aa = mnk_step1(qw, PM)

mean(abs2,PM.-yy1)


function mnk_step1(xx, yy)
    nn = size(xx,1)
    xx = vcat(xx,ones(1,size(xx,2)))
    AA = xx*xx'
    aa = zeros(nn+1,nn)
    for i=1:nn
        BB = xx*yy[i,:]
        aa[:,i] = AA\BB
    end
    yyr = aa'*xx
    return yyr, aa
end


function mnk_step2(xx, yy1, yy, aa)
    nn = size(yy1,1)
    xx = vcat(xx,ones(1,size(xx,2)))
    bb = zeros(nn-1,nn)
    yyr = similar(yy)
    for i = 1:nn
        PT = copy(yy1)
        PT = PT[1:end .!=i,:]
        #PT[i,:] .= 1.0
        AA = PT*PT'
        BB = PT*(yy[i,:].-(aa[:,i]'*xx)')
        bb[:,i] = AA\BB

        yyr[i,:] = bb[:,i]'*PT + aa[:,i]'*xx
    end
    return yyr
end

yy2 = mnk_step2(qw, yy1, PM, aa)
plt = lineplot(yy2[iw,:])
    lineplot!(plt,PM[iw,:])
    println(plt)

mean(abs2, PM.-yy1, dims=2)
mean(abs2,PM.-yy2, dims = 2)

V = cov(PM.-yy2, dims = 2)
W = inv(V)

qw1 = vcat(qw,ones(1,size(qw,2)))
bb = (qw*qw'*W)\(PM*qw'*W)
yy3 = bb*qw
mean(abs2,PM.-(yy2.+yy3), dims = 2)
iw = 2
    plt = lineplot(0*yy2[iw,:].+yy3[iw,:])
    lineplot!(plt,PM[iw,:])
    println(plt)

    for iw = 1:9
        plt = lineplot(yy1[iw,:])
        lineplot!(plt,PM[iw,:])
        println(plt)
    end

plot_P_lr(PM, yy1, 9)

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
