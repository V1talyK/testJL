using ForwardDiff
using Flux, Flux.Tracker, UnicodePlots
include("pde_fun.jl")
include("../grid/rbf.jl")

rw = 0.05/1000;    #Радиус скважины
pw = [50,50]/100;
pk = 1
wxy = [[250, 750],
       [750, 250]]/1000

       #Забойное давление
# TD = Tracker.data;
# xy = [[i,j]/20 for i in 0:20, j in 0:20][:] # Сетка
# wxy = [[0.5 0.5],
#       [0.25 0.25]]
#
 pw = [50.]/100;
 wxy = [[0.61, 0.61]]

bnd, xa, xb, ya, yb = makeModParam()

xy = collect(Iterators.product(10:20:1000, 10:20:1000))
xy = convert(Array{Tuple{Float64,Float64},2},xy)
xy = map(x->[x[1],x[2]],xy)[:]
xy = xy/1000

include("pde_flux.jl")

@time train_lap!(50)


get_hes(x->x[1].^2 +x[2].^3,[1.,2.])

hes_out = get_hes.(x->m3(x),xy)[220:222]
get_hes(x->m3(x),xy[1])

vcat(xy'...).-wxy[1]'

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


loss_1() = Flux.mse(map(x->x[1],m11.(xy)), 0*ones(100))

loss_flux()
loss_flux3()

m4(xy[1])
fRBF(xy[1]')
pde_trialA(xy[1], m1(xy[1]))
bnd = [0 1000; 0 1000]
xy1 = map(x->[(x[1]-bnd[1,1])/(bnd[1,2]-bnd[1,1]),(x[2]-bnd[2,1])/(bnd[2,2]-bnd[2,1])],xy)
funB0(xy1)

prod(pw .-R(xy[1],wxy,pw,rw))

xyg = hcat(map(x->[x[1],x[2]],collect(Iterators.product(0:50:1000,0:50:1000))[:])...)'
Z = fRBF(xyg)
plt = heatmap(reshape(Z,21,21), xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);   display(plt);
xyg[83,:]
