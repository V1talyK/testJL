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

# Example: Minimize Rosenbrock function
function rosenbrock(x)
    return (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
end

function ∇rosenbrock(x)
    return [
        -2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
        200 * (x[2] - x[1]^2)
    ]
end

x0 = [0.0000, 0.0]
x_opt = nonlinear_conjugate_gradient(rosenbrock, ∇rosenbrock, x0, β_method="FR")
println("Optimal solution: x = ", x_opt)
println("Minimum value: f(x) = ", rosenbrock(x_opt))

rosenbrock([1, 1])