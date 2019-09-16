fun23(x) = 0.
Ax(x) = x[2] * sin(pi * x[1]);
psy_trial(net_out,x) = Ax(x) .+ x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * net_out
psy_trial(x) = Ax(x) .+ x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * m1(x)

cb1 = ()->println(loss_flux())
cb2 = ()->println(loss_1())

opt = ADAM(0.1);#opt = Descent(0.0001)

data = Iterators.repeated((), 10)

m1 = Chain(Dense(2,10,σ),Dense(10,1))
m3(x) = psy_trial(x)[1]
m1.(xy)
m3.(xy)

df(x) = Tracker.gradient(x->m1(x)[1],x; nest = true)[1]
d2f(x) = Tracker.gradient(df(x),x; nest = true)
df(xy[1])
d2f(xy[1])

m11(x) = m1(x)[1]
Tracker.hessian(m11,xy[1])
# d2f(xy[1])

Tracker.gradient(m1,xy[1]; nest = true)

function loss_flux()
    #hes_out = ForwardDiff.hessian.(x->Tracker.data(m3(x))[1],xy)
    hes_out = m3.(xy)#Tracker.gradient.(m3(x))[1],xy)
    d2P_dx2 = map(x->x[1],hes_out)
    d2P_dy2 = map(x->x[1],hes_out)
    r_part = fun23.(xy)
    l_part = d2P_dx2 + d2P_dy2;
    B = sum(abs2.(l_part-r_part)) # loss function
    Tracker.TrackedReal{Float64}(B)
end

loss_1() = Flux.mse(map(x->x[1],m11.(xy)), 0*ones(100))

loss_flux()
pr = params(m3)
prm = Flux.params(m1)
prm = Flux.params(m3)

Flux.train!(loss_flux, prm, data, opt, cb = cb1)
Flux.train!(loss_1, prm, data, opt, cb = cb2)

grd = Flux.Tracker.gradient(()->loss_flux(), prm);
grd[W2]

ForwardDiff,gradient(loss_flux,W1)
Flux.Tracker.update!(opt, prm, grd)
prm

NeS = reshape(Tracker.data.(m3.(xy)),10,10)
plt = heatmap(AnS, xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);
plt = heatmap(NeS, xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);
plt = heatmap(PyS, xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);
display(plt)

ps = Params(prm)
Tracker.gradient(ps)
gr
