function fun1(Ab,an,j)
    local x
    for i=1:10
        an[:]=Ab[1]\Ab[2];
    end
end

function fun2(Ab,an,j)
    local x
    for i=1:10
        #display(i)
        #APLU=LinearAlgebra.lufact(Ab[1]);
        #A_ldiv_B!(APLU, Ab[2])
        #an[:]=LU\Ab[2];
        #Ab[1]\Ab[2]
        IterativeSolvers.cg(Ab[1], Ab[2])
    end
end

function fun3(Ab,an,j)
    local A1, a1
    A1 = -Ab[1];
    for i=1:10
        #display(i)
        #display(typeof(A))
        #display(issymmetric(-Ab[1]))
        CL=cholfact(A1);
        #an[:]= -(CL\Ab[2])
        #an[:]=LU\Ab[2];
    end
end

A=rand(1000,1000);
b = rand(1000)
Av=Vector(undef,8); for i=1:8 Av[i]=A; end;
bv=Vector(undef,8); for i=1:8 bv[i]=b; end;
Ab=Vector(undef,8); for i=1:8 Ab[i]=[A,b]; end;

@time fun1(A,b)
iN=Vector(undef,8); for i=1:8 iN[i]=i; end;

@time map(fun1,Ab,an,iN)
@time map(fun2,Ab,an,iN)
@time map(fun3,Ab,an,iN)

@time tmap(fun1,Ab,an,iN)
@time tmap(fun2,Ab,an,iN)
@time tmap(fun3,Ab,an,iN)

@time asyncmap(x->Task(fun2(x)),Ab; ntasks=4)
@time pmap(fun2,Ab;distributed=true)
@time pmap(fun2,Ab)

Ab2 = copy(Ab)
an2 = copy(an)

ai = NaN*zeros(Int32,8)
an = Vector(undef,8); for i=1:8 an[i]=zeros(Int32,length(b)); end;
at = Vector(undef,8); for i=1:8 an[i]=zeros(Int32,length(b)); end;
@time @inbounds Threads.@threads for i=1:6
    ai[i] = Threads.threadid()
    #if Threads.threadid()==1
    #sleep(0.001)
    fun2(Ab[i],an[i],i);
    #end
    #if Threads.threadid()==2
        #fun3(Ab2[i],an2[i],i);
    #end
    #println(i)
end

@time for i=1:8
    fun2(Ab[i],an[i],i);
end

function g(a, n)
   @inbounds Threads.@threads for i in 1:n
       a[i] = i*rand() + rand()^2 * (i - 1)^2 - 2;
   end
   return
end

n = 100000; a = zeros(Float64, n); @time g(a, n)


function tmap(f,arg...)
    y=Vector{Any}(length(arg[1]))
    @inbounds Threads.@threads for i=1:length(arg[1])
        y[i] = f(map(x->x[i],arg));
    end
    return y
end
ccall((:openblas_get_num_threads64_, Base.libblas_name), Cint, ())
