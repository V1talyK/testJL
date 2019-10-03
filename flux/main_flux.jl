#using ForwardDiff
using Flux, Flux.Tracker, UnicodePlots
include("pde_fun.jl")
include("pde_flux.jl")

rw = 0.05;    #Радиус скважины
pw = [50,70];
wxy = [[250 250],
       [750 750]]

        #Забойное давление
TD = Tracker.data;
xy = [[i,j]/20 for i in 0:20, j in 0:20][:] # Сетка
wxy = [[0.5 0.5],
       [0.25 0.25]]

pw = [10.];
wxy = [[0.5 0.5]]

xa=[100,-100];
xb = [0.,0.5]
ya = [100,-100]
yb = [0., 0.5]

train_lap!()

println(Axy([0.5,0.5]))
println(b_out([0.5,0.5]))
b_out(x) = x[2] * sin(pi * x[1]);

(x->w0(x)*(1-funBW(x)))([0.5,0.5])


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

dx = 500;
dy = 500;
np = Int64(1000/dx)
xy = vcat(hcat(collect(0:dx:999),fill(0,np)),
        hcat(fill(999,np), collect(0:dy:999)),
        hcat(collect(999:-dx:0),fill(999,np)),
        hcat(fill(0,np), collect(999:-dy:0)))
z = fill(100,size(xy,1))
xy = vcat(xy,[250 250;750 750])
z = vcat(z,pw)

fRBF = interpByRBF(xy,z);

xyg = hcat(map(x->[x[1],x[2]],collect(Iterators.product(0:50:1000,0:50:1000))[:])...)'
Z = fRBF(xyg)
plt = heatmap(reshape(Z,21,21), xscale=0.1, yscale=0.1, xoffset=0, colormap=:inferno);   display(plt);
xyg[83,:]
