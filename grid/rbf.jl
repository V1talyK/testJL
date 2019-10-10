function funRBF(x,w,xy,ε)
    r2 = map(i->ε*sum(d2p(x,xy[i,:])),1:size(xy,1))
    fi = sqrt.(1 .+r2)
    return sum(fi.*w)
end

@inline d2p(x0,x) = (x0.-x).^2

function interpByRBF(xy,z,ε = 1)
    #Возвращает функцию интерполянт обученную на входных данных
    lp = length(z)
    r2 = zeros(lp,lp)
    for i=1:lp
        r2[:,i] = ε*sum(d2p(xy,xy[i,:]'),dims=2)
    end

    ϕ = sqrt.(1 .+r2);
    w = ϕ\z;
    fun(x) = map(i->funRBF(x[i,:],w,xy,ε),1:size(x,1));
    return fun
end

function test_rbf(np=10;plot = false)
    #np - число точек для обучения
    xy = rand(np,2);
    xyP = collect(Iterators.partition(xy', 2))
    z = (x->x[1].^2 + x[2].^2)
    #println(z.(xPy))

    f1 = interpByRBF(xy,z.(xyP));#ε=1

    xy = collect(Iterators.product(0:0.1:1,0:0.1:1))[:]
    xy = map(x->[x[1],x[2]],xy)[:]
    z_grid = z.(xy)
    z_rbf = f1(hcat(xy...)')
    er = sum(abs.(z_grid-z_rbf))
    if plot
        plt = heatmap(reshape(z_rbf,11,11), xscale=0.1, yscale=0.1,
                 xoffset=0, colormap=:inferno);
        display(plt);
    end
    return er
end

test_rbf(10,plot=false)