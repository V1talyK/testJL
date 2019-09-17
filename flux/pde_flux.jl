cb1 = ()->println(loss_flux())
cb2 = ()->println(loss_1())

opt = ADAM(0.01);#opt = Descent(0.0001)

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
    B = sum(abs2.(l_part-r_part)) # loss function
    #Tracker.TrackedReal{Float64}(B)
end

loss_1() = Flux.mse(map(x->x[1],m11.(xy)), 0*ones(100))

loss_flux()
pr = params(m3)
prm = Flux.params(m1)
prm = Flux.params(m3)

Flux.train!(loss_flux, prm, data, opt, cb = cb1)

NeS = (x->Tracker.data(m3(x))[1]).(xy)
NeS = reshape(NeS,21,21)
plt = heatmap(AnS, xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);   display(plt);
plt = heatmap(NeS, xscale=0.01, yscale=0.01, xoffset=0, colormap=:inferno);    display(plt);
plt = heatmap(PyS, xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);     display(plt);

plt = heatmap(clamp.(NeS,0,1).-AnS, xscale=0.01, yscale=0.01, xoffset=0, colormap=:inferno); display(plt)

NeS[11,:]
P = AnS[11,:]
plt = lineplot(1:21, NeS[11,:], title = "P", name = "w", xlabel = "Sw", ylabel = "kri");
lineplot!(plt,1:21,P)
display(plt)
