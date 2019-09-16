using Pkg
Pkg.add("Knet")
using Knet, StatsBase

struct Linear; w; b; end                  # new type that can be used as a function
(f::Linear)(x) = f.w * x .+ f.b
(f::Linear)() = true           # prediction function if one argument
(f::Linear)(x,y) = mean(abs2, f(x) .- 1)


m1 = Linear(Param(rand(10,2)),Param(rand(10)))
m2 = Linear(Param(rand(1,10)),Param(rand(1)))

M1(x) = psy_trial(m2(m1(x)),x);

fun23(x) = 0.
Ax(x) = x[2] * sin(pi * x[1]);

function psy_trial(net_out,x)
    B = Ax(x) .+ x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * net_out
    return B
end
function loss_knet(xy)
    xy = convert.(Array{Float64,1},xy)
    hes_out = ForwardDiff.jacobian.(M1,xy)#ForwardDiff.hessian.(m1,xy)
    d2P_dx2 = map(x->x[1],hes_out)
    d2P_dy2 = map(x->x[2],hes_out)
    r_part = fun23.(xy)
    l_part = d2P_dx2 + d2P_dy2;
    B = sum(abs2.(l_part-r_part)) # loss function
end

loss_knet(xy)
l1 = @diff loss_knet(xy)
g1 = grad(l1, m1.w)
m0(x) = m2(m1(x))
m0(xy[1])
loss(x) = sum(m0.(x)[1].-1);
loss(xy)
l1 = @diff loss(xy);
g1 = grad(l1,m1.w)


ForwardDiff.hessian.(M1,xy)
ForwardDiff.jacobian.(M1,xy)

dM2(x) = @diff M1(x)
@diff dM2(xy[1])
