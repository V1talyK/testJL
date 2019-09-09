using Flux
using Flux.Tracker

abc = param(abc)
function prbl(x,abc)
    a, b, c = abc;
    y = a.*x.^2 .+b*x .+ c
end

function loss()
    y = prbl(x,p)
    sum(abs2.(y0-y)) # loss function
end

x0 = rand(50); y0 = (x.-0.5).^2 # Dummy data

abc = [1,2,3]
loss()



m = Chain(Dense(1, 2, σ),Dense(2, 1, σ))

p = param([0,0,0])
params = Flux.Params([p])

opt = ADAM(0.1);#opt = Descent(0.1)

data = Iterators.repeated((), 100)
for i=1:5
 Flux.train!(loss, params, data, opt)#cb = ()->println(loss(x,y))
end


plt = scatterplot(x, y, title = "ET", name = "a=0.1", xlabel = "x", ylabel = "y");
scatterplot!(plt,x2, y1, name = "m(x)");

#plt = scatterplot(y1, y, title = "ET", name = "a=0.1", xlabel = "x", ylabel = "y");
display(plt)
