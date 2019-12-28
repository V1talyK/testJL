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
for i=1:10
    c[i] = fun2(i)
end

for i=1:10
    a[i] = fun1()
end

function fun1(tsk)
    wait(tsk);
    b = fetch(tsk);
    return b+1
end


function fun2(i)
    sleep(2)
    return ones(5)*i
end
