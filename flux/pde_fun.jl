fun23(x) = 0#(x[1]==0.5) & (x[2].==0.5) ? 1 : 0;
function Axy(xy)
    #b_out = x[2] * sin(pi * x[1]);
    x = xy[1];    y = xy[2];

    s1 = f0(y)
    s2 = f1(y)
    s3 = g0(x)-((1-x)*g0(0)+x*g0(1))
    s4 = g1(x)-((1-x)*g1(0)+x*g1(1))

    b_out = (1-x)*s1+x*s2+(1-y)*s3 + y*s4
    return b_out
end

# function AxyW(xy,wxy,rw=0.05)
#     F = Axy(xy)
#     F0 = Axy.(wxy)
#     Rw = R(xy,rw)
#     OBW = outBnd.(wxy)
#     B = F - (F0-Rw)*outBnd(xy)/OBW)
#     return B
# end

function AxyW(xy,wxy,pw)
    F = Axy(xy)
    F0 = Axy.(wxy)
    #Rw = R(xy,rw)
    OBW = outBnd.(wxy)
    z = outBnd(xy)

    nw = length(wxy)

    MB = OBW.-OBW';
    MB[1:nw+1:nw*nw] .= 1
    BZ = prod(MB,dims=2)[:]

    AZ = z.-OBW';
    AX = ones(nw)
    AZ1 = ones(nw)
    for i=1:nw
        AZ1[:]=AZ[:]
        AZ1[i]=1.
        AX[i] = prod(AZ1)
    end
    #AZ[1:nw+1:nw*nw] .= 1

    DX = z.*AX./OBW./BZ
    Rw = R(xy,wxy,pw,rw)
    B = F -sum((F0 .- Rw).*DX)
    return B
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

f0(x) = 10+log(rw/sqrt((x-0.5)^2+0.25))
f1(x) = 10+log(rw/sqrt((x-0.5)^2+0.25))
g0(x) = 10+log(rw/sqrt((x-0.5)^2+0.25))
g1(x) = 10+log(rw/sqrt((x-0.5)^2+0.25))#1-(x-0.5).^2
w0(x) = 1


outBnd(xy, xy0 = [0 0; 1 1]) = prod(xy'.-xy0)

function pde_trialA(x, NeIn)

    ψ = funB0(x)*prod(pw .-R(x,wxy,pw,rw))*NeIn
    return ψ
end

function pde_trialB(x, NeIn)
    ψ = Bxy(x) .+ funB0(x)/(1-x[2])*(NeIn-NeIn(x,1)-dN_dy(x,1))s
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

function one_well(xy,rw=0.05)s
    xy0 = [0.5, 0.5]
    R = sum((xy.-xy0).^2).^0.5
    #println(R)
    P = R.==0 ? pw[1] : pw[1]+log(rw/R)
    return P
end

function R(x,wxy,pw,rw=0.05)
    B = zeros(length(pw))
    for i=1:length(pw)
        B[i] = pw[i] +log(rw/sqrt(sum((x.-wxy[i]).^2) + rw*rw))
    end
    return B
end


function loss_flux2()
    hes_out = ForwardDiff.hessian.(x->Tracker.data(m3(x))[1],xy)
    d2P_dx2 = map(x->x[1],hes_out)
    d2P_dy2 = map(x->x[4],hes_out)
    r_part = fun23.(xy)
    l_part = d2P_dx2 + d2P_dy2;
    B = sum(abs2.(l_part-r_part)) # loss function
    #Tracker.TrackedReal{Float64}(B)
end

a=100;
b = 0.25;
si(x,a,b) = 1/(1+exp(-(a*(x-b))))
siab(xy,a,b) = (x->si(x[1],a,b)).(xy)
kh = siab(xy,a,0.25).*(1 .-siab(xy,a,0.75))
dkh =
function tuneKH!(kh,xy)

    for (k, v) in enumerate(xy)
        kh[k] = .&(0.25<v[2]<0.75, 0.25<v[1]<0.75) ? 2 : 1
    end
end
