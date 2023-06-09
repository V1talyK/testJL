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

iw = 1
yy2, aa = two_step(qw, PM) #Прокси
yy2, aa = two_step(qw, ppl2) #Полноценная

plt = lineplot(yy2[iw,:])
    lineplot!(plt,ppl2[iw,:])
    println(plt)

plot_P_lr(ppl2, yy2, 9)

mean(abs2, PM.-yy1, dims=2)
mean(abs2, PM.-yy2, dims = 2)

mnk_step3(qw, yy2, ppl2, aa)

V = cov(ppl2.-yy2, dims = 2)
cov(ppl2[1,:].-yy2[1,:])
cov(ppl2[2,:].-yy2[2,:])
e1 = ppl2[1,:].-yy2[1,:]
e2 = ppl2[2,:].-yy2[2,:]
e1'*e2/33
ee = ppl2.-yy2

V = ee*ee'./33
W = inv(V)

X1 = copy(hcat(qw',ones(33)))
X1 = hcat(X1,zeros(33,10*8))
for i = 2:9
    tmp = hcat(zeros(33,10*(i-1)),hcat(qw',ones(33)),zeros(33,10*(9-i)))
    X1 = vcat(X1, tmp)
end

VV = kron(V,diagm(0=>ones(33)))
bsur = (X1'*inv(VV)*X1)\(X1'*inv(VV)*ppl2'[:])

X1*bsur
mape(ppl2'[:],yy2'[:])
mape(ppl2'[:],X1*bsur)

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
