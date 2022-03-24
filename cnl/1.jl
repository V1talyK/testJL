n = 100000
a = rand(n)
b = rand(n)
c = rand(n)
d = rand(n)
x = rand(n)
x1 = rand(n)


vc = 1:n


@time vmuldiffadd!(view(d,vc),view(c,vc),view(a,vc),view(b,vc))
@time vmuldiffadd!(d,c,a,b)
@profiler for i=1:100 vmuldiffadd!(view(d,vc),view(c,vc),view(a,vc),view(b,vc)) end;

@time d.=d.+c.*(a.-b);
@time for i=1:10
    foo!(a,b,c,d)
end

vx, cv, av, bv = view(x1,vc),view(c,vc),view(a,vc),view(b,vc)

@time for i=1:10
    #BLAS.copyto!(x1,view(d,vc))
    #BLAS.blascopy!(n, d, 1, x1, 1)
    foo!(vx,cv,av,bv)
    #foo!(view(x1,vc),view(c,vc),view(a,vc),view(b,vc))
end

function foo!(a,b,c,d)
    #println(1)
    #axpy!(-1,b,a)
    broadcast!(-,b,a,b);
    broadcast!(*,b,c,b);
    broadcast!(+,d,d,b);
end
