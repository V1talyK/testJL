cb1 = ()->println(loss_flux3())
cb2 = ()->println(loss_1())

m1 = Chain(Dense(2,10,σ),Dense(10,1))
function M3(xy,wxy,pw)

    fun = function m3(x)
        A  = AxyW(x,wxy,pw)
        #ia = findall((x[1].==map(x->x[1],xy)) .& (x[2].==map(x->x[2],xy)))
        A.+pde_trialA(x,m1(x))[1]
    end
    return fun
end

m4 = M3(xy,wxy,pw)

ib = map(x-> findall(sum((vcat(xy'...).-x).^2,dims=2)[:].<rw^2)[1],wxy)
funK, dk_dx, dk_dy = funKH(xa,xb,ya,yb);
funKH(x) = funK(x[1],x[2])*9. +1.;
f_dk_dx(x) = dk_dx(x[1],x[2])*9
f_dk_dy(x) = dk_dy(x[1],x[2])*9

function common_loss(ib)
    kh = funKH.(xy)
    dkdx = f_dk_dx.(xy)
    dkdy = f_dk_dy.(xy)
    lf = function loss_flux()
        #hes_out = ForwardDiff.hessian.(x->Tracker.data(m3(x))[1],xy)
        #hes_out = ForwardDiff.hessian.(x->m3(x,TD(m1(x)))[1],xy)
        #hes_out = m3.(xy)#Tracker.gradient.(m3(x))[1],xy)
        out= get_hes.(x->m4(x),xy)
        #grad_out, hes_out
        dP_dx = map(x->x[1][1],out)
        dP_dy = map(x->x[1][2],out)

        d2P_dx2 = map(x->x[2][1],out)
        d2P_dy2 = map(x->x[2][2],out)
        r_part = fun23.(xy)
        l_part = kh .* (d2P_dx2 .+ d2P_dy2) + (dkdx.*dP_dx + dkdy.*dP_dx);
        l_part[ib].=0
        B = sum(abs2.(l_part-r_part)) # loss function
        #Tracker.TrackedReal{Float64}(B)
    end
    return lf
end

loss_flux3 = common_loss(ib)

function train_lap!()
    opt = ADAM(0.2);#opt = Descent(0.0001)
    data = Iterators.repeated((), 50)

    prm = Flux.params(m1)
    Flux.train!(loss_flux3, prm, data, opt, cb = cb1)
end
