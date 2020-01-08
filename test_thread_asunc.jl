using LinearAlgebra
t1 = Threads.@spawn for i=1:100 lu(A)\b; end;
t2 = Threads.@spawn for i=1:100 lu(A1)\b1; end;
fetch(t1)
fetch(t2)

Threads.nthreads()

A\b
A1 = deepcopy(A)
b1 = deepcopy(b)
Base.Threads.SpinLock()

using Pkg
Pkg.clone("https://github.com/JuliaSmoothOptimizers/HSL.jl.git")

c = Vector(undef,10);
a = Vector(undef,10);
@time for i=1:10
    c[i] = @task fun1(i)
end

a0 = zeros(5);
for i=1:10
    a[i] = fun2(a0,c[i],i)
    a0.=a[i]
end

function fun2(a0,tsk,i)
    println(i)
    #wait(tsk);
    schedule(tsk)
    b = fetch(tsk);
    #println(-i)
    return a0.+1 .+b
end


function fun1(i)
    sleep(2)
    return ones(5)*i/10
end



istaskstarted(c[1])
schedule(c[1]);
istaskdone.(c)
fetch.(c)
