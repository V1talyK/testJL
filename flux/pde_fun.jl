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
    OBW = (x->outBnd(x,[0 0; 1000 1000])).(wxy)
    z = outBnd(xy,[0 0; 1000 1000])

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
    df = Vector(undef,nw)
    for i=1:nw
        AZ1[:]=AZ[:]
        AZ1[i]=1.
        AZ2 = OBW[i].-OBW;
        AZ2[i] = 1.
        df[i] = AZ1[:]./AZ2
        df[i][isnan.(df[i])].=1.
        df[i][isinf.(df[i])].=1.
    end
    #AZ[1:nw+1:nw*nw] .= 1
    prdf = prod.(df)
    #DX = z.*AX./OBW./BZ
    DX = z.*prdf
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

f0(x) = 100;#10+log(rw/sqrt((x-0.5)^2+0.25))
f1(x) = 100;#10+log(rw/sqrt((x-0.5)^2+0.25))
g0(x) = 100;#10+log(rw/sqrt((x-0.5)^2+0.25))
g1(x) = 100;#10+log(rw/sqrt((x-0.5)^2+0.25))#1-(x-0.5).^2
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
    A = zeros(eltype(f2),length(x))
    B = zeros(eltype(f2),length(x))
    xmdx = similar(x)
    xpdx = similar(x)
    for (i, v) in enumerate(x)
        xmdx[:] = x; xmdx[i] = x[i] - dx[i]
        xpdx[:] = x; xpdx[i] = x[i] + dx[i]
        A[i] = ((f(xpdx)-f(xmdx))/2/dx[i])[1]
        B[i] = ((f(xmdx)-f2+f(xpdx))./dx[i].^2)[1]
    end
    return A, B
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

function funKH(xa,xb,ya,yb)
    if length(xa)>0
        fsiX = map(y->x->si(x,y[1],y[2]),zip(xa,xb))
        fsipX = (x->prod(map(f->f(x),fsiX)))
        fX = (x->sum(map(f->1. -f(x),fsiX)))
    else
        fsipX(x) =1;
        fX(x) = 0
    end

    if length(ya)>0
        fsiY = map(y->x->si(x,y[1],y[2]),zip(ya,yb))
        fsipY = (y->prod(map(f->f(y),fsiY)))
        fY = (y->sum(map(f->1. -f(y),fsiY)))
    else
        fsipY(y) = 1;
        fY(y) = 0
    end

    funK(x,y)=fsipX(x)*fsipY(y);
    dk_dx(x,y) = funK(x,y)*fX(x)
    dk_dy(x,y) = funK(x,y)*fY(y)

    return funK, dk_dx, dk_dy
end

si(x,a,b) = 1/(1+exp(-(a*(x-b))))
