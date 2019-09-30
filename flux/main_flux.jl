#using ForwardDiff
using Flux, Flux.Tracker, UnicodePlots
include("pde_fun.jl")
include("pde_flux.jl")

rw = 0.0005;    #Радиус скважины
pw = [10,10];         #Забойное давление
TD = Tracker.data;
xy = [[i,j]/20 for i in 0:20, j in 0:20][:] # Сетка
wxy = [[0.5 0.5],
       [0.25 0.25]]

pw = [10.];
wxy = [[0.5 0.5]]

xa=[100,-100];
xb = [0.25,0.75]
ya = []
yb = []

train_lap!()

println(Axy([0.5,0.5]))
println(b_out([0.5,0.5]))
b_out(x) = x[2] * sin(pi * x[1]);

(x->w0(x)*(1-funBW(x)))([0.5,0.5])


get_hes(x->x[1].^2 +x[2].^3,[1.,2.])

hes_out = get_hes.(x->m3(x),xy)[220:222]
get_hes(x->m3(x),xy[1])



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
