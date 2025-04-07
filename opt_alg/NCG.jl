using LinearAlgebra

function nonlinear_conjugate_gradient(f, ∇f, x0; β_method="FR", max_iter=1000, tol=1e-6)
    x = copy(x0)
    g = ∇f(x)  # Current gradient
    p = -g     # Initial search direction
    f_prev = f(x)
    
    for k in 1:max_iter
        # Line search (using backtracking Armijo condition)
        alp= 1.0
        c = 1e-4  # Armijo condition constant
        ro = 0.5   # Backtracking factor
        while f(x + alp* p) > f(x) + c * alp* dot(g, p)
            alp*=ro
            if alp < 1e-10
                break
                #error("Line search failed: step size too small")
            end
        end
        
        # Update x
        x_new = x + alp* p
        g_new = ∇f(x_new)
        
        # Check convergence
        if norm(g_new) < tol
            println("Converged in $k iterations")
            return x_new
        end
        
        # Compute β (Fletcher-Reeves, Polak-Ribière, or Hestenes-Stiefel)
        if β_method == "FR"
            β = dot(g_new, g_new) / dot(g, g)
        elseif β_method == "PR"
            β = max(0, dot(g_new, g_new - g) / dot(g, g))
        elseif β_method == "HS"
            β = dot(g_new, g_new - g) / dot(p, g_new - g)
        else
            error("Unknown β update method: choose 'FR', 'PR', or 'HS'")
        end
        
        # Update search direction
        p = -g_new + β * p
        
        # Update for next iteration
        x = x_new
        g = g_new
        f_prev = f(x)
    end
    
    println("Maximum iterations reached")
    return x
end



function bounded_ncg(f, ∇f, x0, lb, ub; β_method="FR",
                    max_iter=1000,
                    tol=1.0f-6,
                    xtol=1.0f-6,
                    alp_init=1.0f0,
                    c=1.0f-4,
                    ro=0.5f0)
    
    x = clamp.(x0, lb, ub)
    g = zeros(eltype(x0), length(x))
    p = zeros(eltype(x0), length(x)) # Initialize search direction
    bounded_ncg!(x, g, p, f, ∇f, x0, lb, ub; β_method="FR",
        max_iter=max_iter,
        tol=tol,
        xtol=xtol,
        alp_init=alp_init,
        c=c,
        ro=ro)
    return x
end

function bounded_ncg!(x, g, p, f, ∇f0, x0, lb, ub; β_method="FR", 
                                        max_iter=1000,
                                        tol=1f-6,
                                        xtol=1.0f-6,
                                        alp_init=1f0,
                                        c=1f-4,
                                        ro=0.5f0)
    # Project initial guess to feasible region

    f32max = sqrt(floatmax(Float32))/10;
    ∇f(x)  = clamp.(∇f0(x),-f32max,f32max)
    g_new = copy(g)
    x_new = copy(g)
    g1 = copy(g)
    g .= ∇f(x)

   # Compute projected gradient for initial point
    projected_gradient!(g, x, lb, ub)
    p .= .-g  # Initial search direction

    for k in 1:max_iter
        # Perform projected line search
        alp = projected_backtracking(x_new, g1, g, f, ∇f, x, p, lb, ub; alp_max=alp_init, c=c, ro=ro)
        #alp = two_way_backtracking(x_new, pg, g, f, ∇f, x, p; alp_init=alp_init, ro=ro)
        # Update position with projection
        x_new .= clamp.(x .+ alp .* p, lb, ub)
        if maximum(abs, x_new.-x)<xtol
            #println("xtol reached in $k iterations")
            break
        end

        g_new .= ∇f(x_new)
        projected_gradient!(g_new, x_new, lb, ub)
        # Check convergence
        if norm(g_new) < tol
            println("Converged in $k iterations")
            x .= x_new
            return nothing
        end
        
        # Compute β (Fletcher-Reeves, Polak-Ribière, or Hestenes-Stiefel)
        #println("β=",pg_new, " ", pg)
        if β_method == "FR"
            β = clamp(dot(g_new, g_new), -f32max, f32max) / max(1f-10, clamp(dot(g, g), -f32max, f32max))
        elseif β_method == "PR"
            β = max(0, dot(g_new, g_new .- g) / max(1f-10, dot(g, g)))
        elseif β_method == "HS"
            β = dot(g_new, g_new .- g) / max(1f-10, dot(p, g_new .- g))
        else
            error("Unknown β update method: choose 'FR', 'PR', or 'HS'")
        end
        
        # Update search direction
        #println("u=", pg_new," ",β," ",p) 
        p .= β .* p .- g_new
        #p[isnan.(p)] .= 0f0
        #p .= clamp.(p, -f32max, f32max)
        # Update for next iteration
        x .= x_new
        g .= g_new
    end
    
    #println("Maximum iterations reached")
    return nothing
end

function projected_gradient!(g, x, lb, ub)
    for i in eachindex(x)
        if x[i] <= lb[i] && g[i] > 0
            g[i] = 0.0
        elseif x[i] >= ub[i] && g[i] < 0
            g[i] = 0.0
        end
    end
    return nothing
end

function projected_backtracking(x_new, g1, g, f, ∇f, x, p, lb, ub; alp_max=1f0, c=1f-4, ro=0.5f0)
    alp = alp_max
    g1.=∇f(x)
    projected_gradient!(g1, x, lb, ub)
    
    # Compute directional derivative
    d = dot(g1, p)
    fx = f(x);
    # Find feasible alp
    while true
        x_new .= x .+ alp .* p
        x_new .= clamp.(x_new, lb, ub)
        fx_ = f(x_new)
        if fx_ <= fx + c * alp * d
            break
        end
        alp *= ro
        if alp < 1e-10
            break
            #error("Line search failed: step size too small")
        end
    end
    
    return alp
end

function two_way_backtracking(x_new, pg, g, f, ∇f, x, p; alp_init=1f0, c1=1f-4, c2=0.9f0, ro=0.5f0, max_iter=1000)
    """
    Two-way backtracking line search satisfying strong Wolfe conditions.
    
    Parameters:
    - f: Objective function
    - ∇f: Gradient function
    - x: Current point
    - p: Search direction
    - alp_init: Initial step size guess
    - c1: Armijo condition constant (typically 1e-4)
    - c2: Curvature condition constant (typically 0.9 for CG, 0.1 for Newton)
    - ro: Backtracking factor
    - max_iter: Maximum iterations
    
    Returns:
    - alp: Selected step size
    """
    alp = alp_init
    g.=∇f(x)
    projected_gradient!(pg, x, g, lb, ub)
    d = dot(pg, p)  

    function fa(alp::Float32) 
        x_new .= x .+ alp.*p
        f(x_new)
    end
    function dfa(alp::Float32) 
        x_new .= x .+ alp.*p
        g .= ∇f(x_new)
        projected_gradient!(pg, x, g, lb, ub)
        return nothing
    end
    fx = f(x);

    #First try forward (increasing alp)
    for i in 1:max_iter
        ϕ_alp = fa(alp)
        if ϕ_alp > fx + c1*alp*d || (i > 1 && ϕ_alp >= fa(alp/ro))
            break
        end
        alp = alp/ro  # Increase step size
    end
    
    # Then try backward (decreasing alp) if forward didn't satisfy curvature
    for i in 1:max_iter
        ϕ_alp = fa(alp)
        dfa(alp)
        #dϕ_alp = dot(pg, p)
        
        if ϕ_alp <= fx + c1*alp*d && abs(dϕ_alp) <= c2*abs(d)
            return alp  # Strong Wolfe conditions satisfied
        end
        
        alp = alp*ro  # Decrease step size
        if alp < 1e-10
            #@warn "Step size too small, line search failed"
            return alp
        end
    end
    
    return alp
end

# Example: Minimize Rosenbrock function
# function rosenbrock(x)
#     return (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
# end

# function ∇rosenbrock(x)
#     return [
#         -2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
#         200 * (x[2] - x[1]^2)
#     ]
# end

# x0 = [0.0000, 0.0]
# x_opt = nonlinear_conjugate_gradient(rosenbrock, ∇rosenbrock, x0, β_method="FR")
# println("Optimal solution: x = ", x_opt)
# println("Minimum value: f(x) = ", rosenbrock(x_opt))

# rosenbrock([1, 1])

# Bounds: x1 ∈ [-0.5, 1.5], x2 ∈ [-1.0, 2.0]
# x0 = [0.0, 0.0]
# lb = [-0.5, -1.0]
# ub = [0.5, 2.0]

# x_opt = bounded_ncg(rosenbrock, ∇rosenbrock, x0, lb, ub, β_method="PR")
# println("Optimal solution: x = ", x_opt)
# println("Minimum value: f(x) = ", rosenbrock(x_opt))
# println("Gradient at solution: ∇f(x) = ", ∇rosenbrock(x_opt))