mutable struct ADAM
    η::Float64     # Learning rate
    β1::Float64    # Exponential decay rate for 1st moment estimates
    β2::Float64    # Exponential decay rate for 2nd moment estimates
    ϵ::Float64     # Small constant for numerical stability
    m::Vector{Float64}  # 1st moment vector
    v::Vector{Float64}  # 2nd moment vector
    t::Int         # Time step
    
    function ADAM(η=0.001, β1=0.9, β2=0.999, ϵ=1e-8, n_params::Int=0)
        m = zeros(n_params)
        v = zeros(n_params)
        new(η, β1, β2, ϵ, m, v, 0)
    end
end

function update!(optimizer::ADAM, params::Vector{Float64}, grads::Vector{Float64})
    optimizer.t += 1
    
    # Update biased first moment estimate
    optimizer.m .= optimizer.β1 .* optimizer.m .+ (1 - optimizer.β1) .* grads
    
    # Update biased second raw moment estimate
    optimizer.v .= optimizer.β2 .* optimizer.v .+ (1 - optimizer.β2) .* grads.^2
    
    # Compute bias-corrected first moment estimate
    m̂ = optimizer.m ./ (1 - optimizer.β1^optimizer.t)
    
    # Compute bias-corrected second raw moment estimate
    v̂ = optimizer.v ./ (1 - optimizer.β2^optimizer.t)
    
    # Update parameters
    params .-= optimizer.η .* m̂ ./ (sqrt.(v̂) .+ optimizer.ϵ)
    
    return params
end

# Initialize optimizer for a model with 10 parameters
optimizer = ADAM(0.001, 0.9, 0.999, 1e-8, 10)

# Dummy parameters and gradients
params = rand(10)
grads = rand(10)

nx = 10
xx = rand(1000,nx)
af = [2.0, -4, 3.1, 5, 1, 3, 5, 1, 2,3]
yy = xx*af

f(x) = sum(abs2, yy .-xx*x)
df(x) = -2.0*xx'*(yy .- xx*x)

    # Initial guess and bounds
p0 = ones(nx).+1
# Update parameters
J = zeros(0)
for i=1:10000
    grads = df(p0)
    p0 .= update!(optimizer, p0, grads)
    #println(p0)
    push!(J, sum(abs2,p0-af))
end

using Plots
plot(J)