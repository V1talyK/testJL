fun23(x) = (x[1]==0.5) & (x[2].==0.5) ? 1 : 0;
function Axy(xy)
    #b_out = x[2] * sin(pi * x[1]);
    x = xy[1];    y = xy[2];

    s1 = f0(y)
    s2 = f1(y)
    s3 = g0(x)-((1-x)*g0(0)+x*g0(1))
    s4 = g1(x)-((1-x)*g1(0)+x*g1(1))
    #s3 = g0(x)
    #s4 = g1(x)
    s5 = w0(x)*(1-funBW(x))*x*y*(1-x)*(1-y)/0.5^4

    b_out = (1-x)*s1+x*s2+(1-y)*s3 + y*s4 + s5
    return b_out
end

function Bxy(x)
    x = xy[1];     y = xy[2];

    s1 = (1-x)*f0(y)
    s2 = x*f1(y)
    s3 = g0(x) - ((1-x)*g0(0) + x * g0(1))
    s4 = g1(x)-((1-x)*g1(0)+x*g1(1))

    b_out = s1+s2+s3 + y*s4
    return b_out
end

f0(x) = log(0.05/sqrt((x-0.5)^2+0.25))
f1(x) = log(0.05/sqrt((x-0.5)^2+0.25))
g0(x) = log(0.05/sqrt((x-0.5)^2+0.25))
g1(x) = log(0.05/sqrt((x-0.5)^2+0.25))#1-(x-0.5).^2
w0(x) = 1

f0(1)
println(Axy([0.5,0.5]))
println(b_out([0.5,0.5]))
b_out(x) = x[2] * sin(pi * x[1]);

(x->w0(x)*(1-funBW(x)))([0.5,0.5])


function pde_trialA(x, NeIn)
    ψ = Axy(x) .+ funB0(x)*funBW(x)*NeIn
    return ψ
end

function pde_trialB(x, NeIn)
    ψ = Bxy(x) .+ funB0(x)/(1-x[2])*(NeIn-NeIn(x,1)-dN_dy(x,1))
    return ψ
end

@inline funB0(x) = prod(x.*(1 .-x))
@inline funBW(x) = sum(1 .-exp.(-(x.-0.5).^2))/2

psy_trial(net_out,x) = Ax(x) .+ x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * net_out
psy_trial(x) = Ax(x) .+ x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * m1(x)

function get_hes(f,x)
    dx = 0.001*x
    dx[x.==0] .= 0.001
    f2 = 2*f(x);
    B = zeros(eltype(f2),length(x))
    xmdx = similar(x)
    xpdx = similar(x)
    for (i, v) in enumerate(x)
        xmdx[:] = x; xmdx[i] = x[i] - dx[i]
        xpdx[:] = x; xpdx[i] = x[i] + dx[i]
        B[i] = ((f(xmdx)-f2+f(xpdx))./dx[i].^2)[1]
    end
    return B
end

get_hes(x->x[1].^2 +x[2].^3,[1.,2.])

hes_out = get_hes.(x->m3(x),xy)
get_hes(x->m3(x),xy[1])

function one_well(xy)
    xy0 = [0.5, 0.5]
    R = sum((xy.-xy0).^2).^0.5
    rw = 0.05
    P = R.==0 ? 1 : log(rw/R)
    return P
end

AnS = one_well.(xy);
AnS = reshape(AnS,21,21)
