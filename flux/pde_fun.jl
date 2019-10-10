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
    F = Axy(xy[613])
    F0 = Axy.(wxy)
    #Rw = R(xy,rw)
    OBW = (x->outBnd(x,[0 0; 1000 1000])).(wxy)
    z = outBnd(xy[613],[0 0; 1000 1000])

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
        AZ2 = Float64.(OBW[i].-OBW);
        AZ2[i] = 1.
        AZ2[AZ2.==0] .= 1e-6
        AZ1[AZ1.==0] .= 1e-6
        df[i] = AZ1[:]./AZ2
        df[i][isnan.(df[i])].=1.
        df[i][isinf.(df[i])].=1.
    end
    #AZ[1:nw+1:nw*nw] .= 1
    prdf = prod.(df)
    #DX = z.*AX./OBW./BZ
    DX = z.*prdf./OBW
    Rw = R(xy[613],wxy,pw,rw)
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
    bnd = [0 1000; 0 1000]
    #x1 = [(x[1]-bnd[1,1])/(bnd[1,2]-bnd[1,1]),(x[2]-bnd[2,1])/(bnd[2,2]-bnd[2,1])]
    ψ = funB0(x)*prod(1 .-fromWell(x,wxy,rw))*NeIn
    #ψ = funB0(x)*prod(fromWell2(x,wxy))*NeIn
    #ψ = fun_GU(x)*prod(pw .-R(x,wxy,pw,rw))*NeIn
    #ψ = funB0(x)*(fRBF0(x')[1])*NeIn

    return ψ
end

function pde_trialB(x, NeIn)
    ψ = Bxy(x) .+ funB0(x)/(1-x[2])*(NeIn-NeIn(x,1)-dN_dy(x,1))s
    return ψ
end

@inline funB0(x,bnd = [0 1;0 1]) = prod((view(bnd,:,1).-x).*(view(bnd,:,2) .-x))
@inline funBW(x) = sum(1 .-exp.(-(x.-0.5).^2))/2

psy_trial(net_out,x) = Ax(x) .+ x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * net_out
psy_trial(x) = Ax(x) .+ x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * m1(x)

function get_hes(f,x)
    dx = 0.001*x;#*ones(length(x))#
    dx[x.==0] .= 0.001
    f2 = 2*f(x);
    A = zeros(eltype(f2),length(x))
    B = zeros(eltype(f2),length(x))
    xmdx = similar(x)
    xpdx = similar(x)
    @inbounds for (i, v) in enumerate(x)
        xmdx[:] = x; xmdx[i] = v - dx[i]
        xpdx[:] = x; xpdx[i] = v + dx[i]
        fp = f(xpdx);
        fm = f(xmdx);
        A[i] = ((fp-fm)/2/dx[i])[1]
        B[i] = ((fm-f2+fp)./dx[i].^2)[1]
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

function fromWell(x,wxy,rw=0.05)
    B = zeros(length(wxy))
    for i=1:length(wxy)
        B[i] = log(rw/(sqrt(sum((x.-wxy[i]).^2)) + rw))
    end
    return B
end

function fromWell2(x,wxy)
    B = sum(d2p(hcat(wxy...),x),dims=1)
    return sqrt.(B')
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

function funKH(xa,xb,ya,yb,fun=si)
    if length(xa)>0
        fsiX = map(y->x->fun(x,y[1],y[2]),zip(xa,xb))
        fsipX = (x->prod(map(f->f(x),fsiX)))
        fX = (x->sum(map(f->f[2]*(1. -f[1](x)),zip(fsiX,xa))))
    else
        fsipX(x) =1;
        fX(x) = 0
    end

    if length(ya)>0
        fsiY = map(y->x->fun(x,y[1],y[2]),zip(ya,yb))
        fsipY = (y->prod(map(f->f(y),fsiY)))
        fY = (y->sum(map(f->f[2]*(1. -f[1](y)),zip(fsiY,ya))))
    else
        fsipY(y) = 1;
        fY(y) = 0
    end

    funK(x,y)=fsipX(x)*fsipY(y);
    dk_dx(x,y) = funK(x,y)*fX(x)
    dk_dy(x,y) = funK(x,y)*fY(y)

    return funK, dk_dx, dk_dy
end

@inline si(x,a,b) = 1/(1+exp(-(a*(x-b))))
@inline isru(x,a,b) = x/sqrt(1+(a*(x-b))^2)

function makeRBFfromBoundary(pk, pw, bnd,wxy, np=10)
    dx = bnd[1][2]/np;
    dy = bnd[2][2]/np;
    maxX, maxY = bnd[1][2]*0.99, bnd[2][2]*0.99;
    xy = vcat(hcat(collect(0:dx:maxX),fill(0,np)),
            hcat(fill(maxX,np), collect(0:dy:maxY)),
            hcat(collect(maxX:-dx:0),fill(maxY,np)),
            hcat(fill(0,np), collect(maxY:-dy:0)))
    z = fill(pk,size(xy,1))
    xy = vcat(xy,hcat(wxy...)')
    z = vcat(z,pw)
    fun = interpByRBF(xy,z,1000);
    return fun
end
