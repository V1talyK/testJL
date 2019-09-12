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

opt = ADAM(0.1);#opt = Descent(0.1)

data = Iterators.repeated((), 10)
W1p = Flux.param(W1)

m11 = Chain(Dense(2,10,σ),Dense(10,1))
m3(x) = psy_trial(m11(x),x)[1]
m3.(xy)

df(x) = Tracker.gradient(m3,x; nest = true)[1][2]
d2f(x) = Tracker.gradient(df,x; nest = false)[1]
# df.(xy)
# d2f.(xy)
# d2f(xy[1])

function loss_flux()
    hes_out = ForwardDiff.hessian.(x->Tracker.data(m3(x)),xy)
    d2P_dx2 = map(x->x[1],hes_out)
    d2P_dy2 = map(x->x[4],hes_out)
    r_part = fun23.(xy)
    l_part = d2P_dx2 + d2P_dy2;
    B = sum(abs2.(l_part-r_part)) # loss function
    Tracker.TrackedReal{Float64}(B)
end
loss_flux()
pr = params(m11)
prm = Flux.params(m11)

Flux.train!(loss_flux, prm, data, opt, cb = cb1)m
cb1 = ()->println(loss_flux())
