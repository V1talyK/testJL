using AutoEncoderToolkit
using Flux

enc = AutoEncoderToolkit.Encoder(Flux.Chain(Dense(784, 400, relu), Dense(400, 20)))
x = rand(Float32, 784);
y = enc(x);
dec = AutoEncoderToolkit.Decoder(Flux.Chain(Dense(20, 400, relu), Dense(400, 784)))

z = dec(y)


η = 1e-3
# Explicit setup of optimizer
opt_1 = Flux.Train.setup(
    Flux.Optimisers.Adam(η),
    vae
)

sum(abs2, x.-z)