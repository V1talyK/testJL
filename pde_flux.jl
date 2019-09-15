function fun23(x)
    return 0.
end

function Ax(x)
    return x[2] * sin(pi * x[1])
end

function psy_trial(net_out,x)
    #x = in_z[1]
    # = in_z[2]
    B = Ax(x) .+ x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * net_out
    return B
end
cb1 = ()->println(loss_flux())
cb2 = ()->println(loss_1())

opt = ADAM(0.1);#opt = Descent(0.0001)

data = Iterators.repeated((), 10)


m11 = Chain(Dense(2,10,σ),Dense(10,1))
m3(x) = psy_trial(m11(x),x)
m1.(xy)

m3.(xy)

df(x) = Tracker.gradient(m3,x; nest = true)[1][2]
d2f(x) = Tracker.gradient(df,x; nest = false)[1]
#df.(xy)
#d2f.(xy)
# d2f(xy[1])

function loss_flux()
    #m3(x) = psy_trial(m11(x),x)[1]
    #hes_out = ForwardDiff.hessian.(x->Tracker.data(m3(x)[1]),xy)
    hes_out = ForwardDiff.hessian.(x->Tracker.data(m1(x)),xy)
    #m23 = x->m11(x)
    #hes_out = ForwardDiff.hessian.(x->Tracker.data(m23(x)[1]),xy)
    d2P_dx2 = map(x->x[1],hes_out)
    d2P_dy2 = map(x->x[4],hes_out)
    r_part = fun23.(xy)
    l_part = d2P_dx2 + d2P_dy2;
    B = sum(abs2.(l_part-r_part)) # loss function
    Tracker.TrackedReal{Float64}(B)
end

loss_1() = Flux.mse(map(x->x[1],m11.(xy)), 0*ones(100))

loss_flux()
pr = params(m3)
prm = Flux.params(m11)
W1,b1,W2,b2 = map(x->param(x),[W1,b1,W2,b2])
prm = Flux.params(W1,b1,W2,b2)

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
