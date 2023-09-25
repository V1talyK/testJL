using Distributions, LsqFit
"Проверяем катомную регресию с пакетом + проверка ошибок"

xx = rand(0:0.1:100,100, 2)
a1 = 5 .* rand(0.95:0.01:1.05,size(xx,1))
a2 = 3 .* rand(0.95:0.01:1.05,size(xx,1))
b1 = -6 .* rand(0.95:0.01:1.05,size(xx,1))

yy = a1 .*xx[:,1] + a2 .*xx[:,2] .-b1
yyr, aa, er_aa = mnk_step1(xx', yy')


sig = sum(abs2,yyr' .- yy)/(length(xx) - length(aa))

CovM3x = inv(Matrix(er_aa))*sig
tp = quantile(TDist(length(xx) - length(aa)),0.975)

Δa = sqrt.(diag(CovM3x))*tp


model(x, p) = p[1] .* x[:,1] .+ p[2] .* x[:,2] .+ p[3]
p0 = [1. ,1., 1.]

fit1 = curve_fit(model, xx, yy, p0);
fit1.param

se = estimate_covar(fit1)
cio = confidence_interval(fit1, 0.05)

map(x->(x[2]-x[1])/2, cio)
