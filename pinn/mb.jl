using NeuralPDE, OrdinaryDiffEq, Lux, Random, OptimizationOptimJL, LineSearches, Plots
using LinearAlgebra

nt = 100
ppl = zeros(Float32, 100)
init_ppl  = 10f0;
pa  = 10f0;
ua = [0.02f0, 0.1f0, 0.1f0]
qw = 10f0*ones(Float32, nt, 1); qw[50:end] .= 5f0;
pc = init_ppl*ones(Float32, nt, 1) 
pc[:] .= collect(range(init_ppl, init_ppl-1, 100))#

mb2!(ppl, ua, qw, pc, pa)
plot(ppl)


function mb(p, u, t)
    p₁ = p[1]
    β₁, a1, a2, pa = u
    dp₁ = -β₁ * q2(t) + a1 * (pcf(t) - p₁) + a2 * (pa - p₁)
    #dp₂ = a12 * (p₁ - p₂) + a2c * (pc - p₂)
    [dp₁]
end

function q2(t)
    return ifelse(t>50, 5.0, 10.0)
end

function pcf(t)
    _, ia = findmin(abs, (1:nt) .- t)
    return pc[ia]
end

tspan = (0.0, 100.0)
p0 = [init_ppl]
init_p = [1.0, 1.0, 1.0, init_ppl]
true_p = vcat(ua, init_ppl)
prob = ODEProblem(mb, p0, tspan, init_p)

prob_data = remake(prob, p = true_p)
sol_data = solve(prob_data, Tsit5(), saveat = 0.1)
t_ = sol_data.t
u_ = reduce(hcat, sol_data.u)

rng = Random.default_rng()
Random.seed!(rng, 0)
n = 15
chain = Chain(Dense(1, n, σ), Dense(n, n, σ), Dense(n, n, σ), Dense(n, 2))
ps, st = Lux.setup(rng, chain) |> f64

additional_loss(phi, θ) = sum(abs2, phi(t_, θ) .- u_) / size(u_, 2)

opt = LBFGS(linesearch = BackTracking())
alg = NNODE(chain, opt, ps; strategy = WeightedIntervalTraining([0.7, 0.2, 0.1], 500),
    param_estim = true, additional_loss)

sol = solve(prob, alg, verbose = true, abstol = 1e-8, maxiters = 50, saveat = t_)

plt = plot(sol_data.t, vcat(sol_data.u...), labels = ["u1_data"] )
    scatter!(plt, ppl)
    plot!(plt, sol, labels = ["u1_pinn"])


sol_data = solve(prob_data, Tsit5(), saveat = 1.0)

ppl1 = zeros(nt)
ppl0 = init_ppl
for t = 1:nt
    ppl1[t] = ppl0 - ua[1]*qw[t] + ua[2]*(pc[t] - ppl0) + ua[3]*(pa - ppl0)
    ppl0 = ppl1[t]
end

plot!(plt, ppl1)