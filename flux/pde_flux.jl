cb1 = ()->println(loss_flux())
cb2 = ()->println(loss_1())

opt = ADAM(0.2);#opt = Descent(0.0001)

data = Iterators.repeated((), 50)

m1 = Chain(Dense(2,10,σ),Dense(10,1))
m3(x) = pde_trialA(x,m1(x))[1]
m1.(xy)
m3.(xy)

df(x) = Tracker.gradient(x->m1(x)[1],x; nest = true)[1]
d2f(x) = Tracker.gradient_nested(x->df(x)[1],x)
df(xy[2])
d2f(xy[2])

m11(x) = m1(x)[1]
Tracker.hessian(m11,xy[1])
# d2f(xy[1])

Tracker.gradient(m1,xy[1]; nest = true)
TD = Tracker.data;
function loss_flux()
    #hes_out = ForwardDiff.hessian.(x->Tracker.data(m3(x))[1],xy)
    #hes_out = ForwardDiff.hessian.(x->m3(x,TD(m1(x)))[1],xy)
    #hes_out = m3.(xy)#Tracker.gradient.(m3(x))[1],xy)
    hes_out = get_hes.(x->m3(x),xy)
    d2P_dx2 = map(x->x[1],hes_out)
    d2P_dy2 = map(x->x[2],hes_out)
    r_part = fun23.(xy)
    l_part = d2P_dx2 + d2P_dy2;
    l_part[221]=0
    B = sum(abs2.(l_part-r_part)) # loss function
    #Tracker.TrackedReal{Float64}(B)
end
loss_1() = Flux.mse(map(x->x[1],m11.(xy)), 0*ones(100))

loss_flux()
loss_flux2()

prm = Flux.params(m1)

Flux.train!(loss_flux, prm, data, opt, cb = cb1)
