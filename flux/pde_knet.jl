using Pkg
Pkg.add("Knet")
using Knet

struct Linear; w; b; end                  # new type that can be used as a function
(f::Linear)(x) = f.w * x .+ f.b
(f::Linear)() = true           # prediction function if one argument
(f::Linear)(x,y) = mean(abs2, f(x) - y)


m1 = Linear(rand(10,2),rand(10))
m2 = Linear(rand(1,10),rand(1))

(x->psy_trial(m2(m1(x)),x)).(xy)
