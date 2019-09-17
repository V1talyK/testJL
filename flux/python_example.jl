using ForwardDiff, UnicodePlots
xy = [[i,j]/10 for i in 1:10, j in 1:10][:]
xy = [[i,j]/20 for i in 0:20, j in 0:20][:]
xy = [[i,j]/100 for i in 0:100, j in 0:100][:]
m1(x) = Tracker.data(psy_trial(m2(xy[1]), x)[1])
m2 = Chain(Dense(2,10,σ),Dense(10,1))


W1 = rand(10,2);
W2 = rand(1,10);
b1 = zeros(10,1);
b2 = zeros(1,1);

modl1(x,W1,b1) = sigmoid.(W1*x.+b1);
modl2(x,W2,b2) = W2*x.+b2;
m1(x) = psy_trial(modl2(modl1(x,W1,b1),W2,b2),x)[1]


Uf.(xy)
ForwardDiff.hessian.(Uf,xy)

m1.(xy))

function analytic_solution(xy)
    x = map(x->x[1],xy)
    y = map(x->x[2],xy)
    B =(1 / (exp(pi) - exp(-pi))) * sin.(pi * x) .* (exp.(pi * y) - exp.(-pi * y))
    return B
end

AnS = reshape(analytic_solution(xy),101,101)

function loss_flux(W1,W2)
    #W1 = W[1]
    #W2 = W[2]
    Uf(x) = psy_trial(modl2(modl1(x,W1,b1),W2,b2),x)[1]
    #net_out = m1.(xy)
    hes_out = ForwardDiff.hessian.(Uf,xy)#ForwardDiff.hessian.(m1,xy)
    d2P_dx2 = map(x->x[1],hes_out)
    d2P_dy2 = map(x->x[4],hes_out)
    r_part = fun23.(xy)
    l_part = d2P_dx2 + d2P_dy2;
    sum(abs2.(l_part-r_part)) # loss function
end

ForwardDiff.hessian.(Uf,xy)
loss_flux(W1,W2)
@time for i=1:50
    loss_grad1 = ForwardDiff.gradient(x->loss_flux(x,W2), W1)
    loss_grad2 = ForwardDiff.gradient(x->loss_flux(W1,x), W2)
    println("lf: $(loss_flux(W1,W2)), as: $(sum(Uf.(xy)[:] .- AnS[:]))")

    W1[:] = W1[:] .- 0.0001*loss_grad1[:]
    W2[:] = W2[:] .- 0.0001*loss_grad2[:]
end


PyS = reshape(Uf.(xy),10,10);
