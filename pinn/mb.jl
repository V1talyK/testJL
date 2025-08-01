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

mb2!(ppl, ua, qw, hcat(pc, fill(pa, nt)), pa)
plt = plot(ppl)


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
init_p = [1f0, 1f0, 1f0, init_ppl]
true_p = vcat(ua, init_ppl)
prob = ODEProblem(mb, p0, tspan, true_p)

prob_data = remake(prob, p = init_p)
sol_data = solve(prob_data, Tsit5(), saveat = 0.1)
ti = collect(1:5:nt)
pf = zeros(Float32, 1, length(ti)); pf.=deepcopy(ppl[ti]')

rng = Random.default_rng()
Random.seed!(rng, 0)
n = 15
chain = Chain(Dense(1, n, σ), Dense(n, n, σ), Dense(n, n, σ), Dense(n, 2))
ps, st = Lux.setup(rng, chain) |> f64

additional_loss(phi, θ) = sum(abs2, phi(ti, θ) .- pf) / length(pf)

opt = LBFGS(linesearch = BackTracking())
alg = NNODE(chain, opt, ps; strategy = WeightedIntervalTraining([0.7, 0.2, 0.1], 50),
    param_estim = true, additional_loss)

sol = solve(prob, alg, verbose = true, abstol = 1e-8, maxiters = 200, saveat = ti)

plt = plot(sol_data.t, vcat(sol_data.u...), labels = ["u1_data"] )
plt = scatter(ppl)
    plot!(plt, sol, labels = ["u1_pinn"])
    scatter!(plt, ti, pf[:])


# sol_data = solve(prob_data, Tsit5(), saveat = 1.0)

# ppl1 = zeros(nt)
# ppl0 = init_ppl
# for t = 1:nt
#     ppl1[t] = ppl0 - ua[1]*qw[t] + ua[2]*(pc[t] - ppl0) + ua[3]*(pa - ppl0)
#     ppl0 = ppl1[t]
# end

# mb2!(ppl, ua, qw, pc, pa, 0.5f0)
# plt = plot(ppl1)
# scatt!(plt, ppl, label = "mb2")


# pplc = zeros(Float32, nt)
# DT = get_index_preess_points0(true, Float32.(ppl1), 1:nt, init_ppl)
# XX1 = create_XX0(DT, hcat(pc, fill(pa, nt)), cumsum(qw, dims = 1), init_ppl, Float32.(ppl1))
# aa = reb_ppl_well_by_mb!(pplc, Float32.(ppl1), 1:nt, init_ppl, hcat(pc, fill(pa, nt)), qw, DT, XX1;
#                         drop_Δt=2)

# plot!(plt, pplc, label = "alg", lw = 2)

# mb2!(ppl, Float32.(aa), qw,  hcat(pc, fill(pa, nt)), pa, 0.5f0)