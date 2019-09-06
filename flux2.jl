using Flux
using Flux.Tracker

W = rand(2, 5)
b = rand(2)

predict(x) = W*x .+ b

loss(x, y) = Flux.mse(m(x), y)

x = rand(50); y = (x.-0.5).^2 # Dummy data
loss(x, y) # ~

m = Chain(Dense(2, 2, σ),Dense(2, 2, σ))

ps = Flux.params(m)
opt = Descent(0.1)

data = [(x[j:j+1],y[j:j+1]) for j=1:49]
data = Iterators.repeated((x, y), 3)
for i=1:100
    Flux.train!(loss, ps, data, opt)#cb = ()->println(loss(x,y))
end

m(x[1:2])

loss(x,y)

x2 = rand(50)
data = [x2[j:j+1] for j=1:49]
y1 = Tracker.data(m(x2))

plt = scatterplot(x2, y1, title = "ET", name = "a=0.1", xlabel = "x", ylabel = "y");
scatterplot!(plt,x, y, name = "m(x)");

#plt = scatterplot(y1, y, title = "ET", name = "a=0.1", xlabel = "x", ylabel = "y");
display(plt)
