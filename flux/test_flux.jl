using Flux, Flux.Tracker
using Flux.Tracker: update!
using Flux: @epochs

W = rand(2, 5)
b = rand(2)

predict(x) = W*x .+ b

function loss(x, y)
  ŷ = predict(x)
  sum((y .- ŷ).^2)
end

x, y = rand(5), rand(2) # Dummy data
loss(x, y) # ~ 3

W = param(W)
b = param(b)

gs = Tracker.gradient(() -> loss(x, y), params(W, b))

Δ = gs[W]
Δ1 = gs[b]

# Update the parameter and reset the gradient
update!(W, -0.1Δ)
update!(b, -0.1Δ1)

loss(x, y) # ~ 2.5

W1 = param(rand(3, 5))
b1 = param(rand(3))
layer1(x) = W1 * x .+ b1

W2 = param(rand(2, 3))
b2 = param(rand(2))
layer2(x) = W2 * x .+ b2

model(x) = layer2(σ.(layer1(x)))

model(rand(5)) # => 2-element vector

function linear(in, out)
  W = param(randn(out, in))
  b = param(randn(out))
  x -> W * x .+ b
end

linear1 = linear(5, 3) # we can access linear1.W etc
linear2 = linear(3, 2)

model(x) = linear2(σ.(linear1(x)))

model(rand(5)) # => 2-element vector

struct Affine
  W
  b
end

Affine(in::Integer, out::Integer) =
  Affine(param(randn(out, in)), param(randn(out)))

# Overload call, so the object can be used as a function
(m::Affine)(x) = m.W * x .+ m.b

a = Affine(10, 5)

a(rand(10)) # => 5-element vector

θ = Params([W, b])
grads = Tracker.gradient(() -> loss(x, y), θ)

opt = Descent(0.1) # Gradient descent with learning rate 0.1

for p in (W, b)
  update!(opt, p, grads[p])
end

loss(x,y)

Flux.train!(loss, (W, b), [(x,y)], opt, cb = ()->println(loss(x,y)))


m = Chain(
  Dense(784, 32, σ),
  Dense(32, 10), softmax)

loss(x, y) = Flux.mse(m(x), y)
ps = Flux.params(m)

# later
x, y = rand(784), rand(10)
@epochs 3 Flux.train!(loss, ps, [(x,y)], opt, cb = ()->display(loss(x,y)))
