m1 = Chain(Dense(2,10,σ),Dense(10,1))
fRBF = makeRBFfromBoundary(pk, pw, bnd,wxy)
fRBF0 = makeRBFfromBoundary(1, 0*pw, bnd,wxy)

function M3(xy,wxy,pw)
    fun = function m3(x)
        A  = fRBF(x')[1];#AxyW(x,wxy,pw)
        #ia = findall((x[1].==map(x->x[1],xy)) .& (x[2].==map(x->x[2],xy)))
        A.+pde_trialA(x,m1(x))[1]
        #A.+m1(x)[1]
    end
    return fun
end

m3 = M3(xy,wxy,pw)

ib = map(x-> findall(sum((vcat(xy'...).-x').^2,dims=2)[:].<rw^2)[1],wxy)
funK, dk_dx, dk_dy = funKH(xa,xb,ya,yb);
funKH(x) = funK(x[1],x[2]) .+0.1;
f_dk_dx(x) = dk_dx(x[1],x[2])
f_dk_dy(x) = dk_dy(x[1],x[2])

function common_loss(ib)
    kh = funKH.(xy)
    dkdx = f_dk_dx.(xy)
    dkdy = f_dk_dy.(xy)
    lf = function loss_flux()
        #hes_out = ForwardDiff.hessian.(x->Tracker.data(m3(x))[1],xy)
        #hes_out = ForwardDiff.hessian.(x->m3(x,TD(m1(x)))[1],xy)
        #hes_out = m3.(xy)#Tracker.gradient.(m3(x))[1],xy)
        out= get_hes.(x->m3(x),xy)
        #grad_out, hes_out
        dP_dx = map(x->x[1][1],out)
        dP_dy = map(x->x[1][2],out)

        d2P_dx2 = map(x->x[2][1],out)
        d2P_dy2 = map(x->x[2][2],out)
        r_part = fun23.(xy)#kh .*
        l_part = (d2P_dx2 .+ d2P_dy2) + 0*(dkdx.*dP_dx + dkdy.*dP_dx);
        l_part[ib].=0
        B = sum(abs2.(l_part-r_part)) # loss function
        #Tracker.TrackedReal{Float64}(B)
    end
    return lf
end

loss_flux = common_loss(ib)
cb1 = ()->println(loss_flux())

function train_lap!()
    #opt = Descent(0.0001)
    opt = ADAM(0.1);#opt = Descent(0.0001)
    data = Iterators.repeated((), 25)

    prm = Flux.params(m1)
    Flux.train!(loss_flux, prm, data, opt, cb = cb1)
end
