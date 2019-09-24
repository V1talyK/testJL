
cb1 = ()->println(loss_flux())
cb2 = ()->println(loss_1())

m1 = Chain(Dense(2,10,σ),Dense(10,1))
m3(x) = pde_trialA(x,m1(x))[1]

ib = map(x-> findall(sum((vcat(xy'...).-x).^2,dims=2)[:].<rw^2)[1],wxy)

function common_loss(ib)
    kh = ones(Float64,length(xy))
    tuneKH!(kh,xy)
    lf = function loss_flux()
        #hes_out = ForwardDiff.hessian.(x->Tracker.data(m3(x))[1],xy)
        #hes_out = ForwardDiff.hessian.(x->m3(x,TD(m1(x)))[1],xy)
        #hes_out = m3.(xy)#Tracker.gradient.(m3(x))[1],xy)
        hes_out = get_hes.(x->m3(x),xy)
        d2P_dx2 = map(x->x[1],hes_out)
        d2P_dy2 = map(x->x[2],hes_out)
        r_part = fun23.(xy)
        l_part = kh .* d2P_dx2 + kh .* d2P_dy2;
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
