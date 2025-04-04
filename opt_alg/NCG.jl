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


function bounded_ncg(f, ∇f, x0, lb, ub; β_method="FR",
                    max_iter=1000,
                    tol=1.0f-6,
                    alp_init=1.0f0,
                    c=1.0f-4,
                    ro=0.5f0)
    
    x = clamp.(x0, lb, ub)
    g = zeros(eltype(x0), length(x))
    p = zeros(eltype(x0), length(x)) # Initialize search direction
    pg = zeros(eltype(x0), length(x)) # Initialize search direction
    bounded_ncg!(x, g, p, pg, f, ∇f, x0, lb, ub; β_method="FR",
        max_iter=max_iter,
        tol=tol,
        alp_init=alp_init,
        c=c,
        ro=ro)
    return x
end

function bounded_ncg!(x, g, p, pg, f, ∇f0, x0, lb, ub; β_method="FR", 
                                        max_iter=1000,
                                        tol=1f-6,
                                        alp_init=1f0,
                                        c=1f-4,
                                        ro=0.5f0)
    # Project initial guess to feasible region

    f32max = sqrt(floatmax(Float32))/10;
    ∇f(x)  = clamp.(∇f0(x),-f32max,f32max)
    pg_new = copy(pg)
    g_new = copy(g)
    pg1 = copy(pg)
    g1 = copy(g)
    g .= ∇f(x)

   # Compute projected gradient for initial point
    projected_gradient!(pg, x, g, lb, ub)
    p .= .-pg  # Initial search direction
    # println("p=",p) 
    # println("g=",g) 
    # println("x=",x) 
    # println("lb=",lb) 
    # println("ub=",ub) 

    for k in 1:max_iter
        # Perform projected line search
        alp = projected_backtracking(pg1, g1, f, ∇f, x, p, lb, ub; alp_max=alp_init, c=c, ro=ro)
        
        # Update position with projection
        #println("x_new=", x," ",alp," ",p) 
        x_new = clamp.(x .+ alp .* p, lb, ub)
        #println("x_new=", x_new) 
        #println(typeof.([x,alp,p]))
        #println("p=",p)
        g_new .= ∇f(x_new)
        #println("x_new=",x_new) 
        #println("g_new=",g_new) 
        #println("g_new=",g_new) 
        #g_new .= clamp.(g_new, -f32max, f32max)
        #println("g_new=",g_new) 
        projected_gradient!(pg_new, x_new, g_new, lb, ub)
        #pg_new .= clamp.(pg_new, -f32max, f32max)
        # println("---------$k")
        # Check convergence
        if norm(pg_new) < tol
            #println("Converged in $k iterations")
            x .= x_new
            return nothing
        end
        
        # Compute β (Fletcher-Reeves, Polak-Ribière, or Hestenes-Stiefel)
        #println("β=",pg_new, " ", pg)
        if β_method == "FR"
            β = clamp(dot(pg_new, pg_new), -f32max, f32max) / max(1f-10, clamp(dot(pg, pg), -f32max, f32max))
        elseif β_method == "PR"
            β = max(0, dot(pg_new, pg_new .- pg) / max(1f-10, dot(pg, pg)))
        elseif β_method == "HS"
            β = dot(pg_new, pg_new .- pg) / max(1f-10, dot(p, pg_new .- pg))
        else
            error("Unknown β update method: choose 'FR', 'PR', or 'HS'")
        end
        
        # Update search direction
        #println("u=", pg_new," ",β," ",p) 
        p .= -pg_new .+ β .* p
        #p[isnan.(p)] .= 0f0
        #p .= clamp.(p, -f32max, f32max)
        # Update for next iteration
        x .= x_new
        g .= g_new
        pg .= pg_new
    end
    
    #println("Maximum iterations reached")
    return nothing
end

function projected_gradient!(pg, x, g, lb, ub)
    for i in eachindex(x)
        if x[i] <= lb[i] && g[i] > 0
            pg[i] = 0.0
        elseif x[i] >= ub[i] && g[i] < 0
            pg[i] = 0.0
        else
            pg[i] = g[i]
        end
    end
    return nothing
end

function projected_backtracking(pg, g, f, ∇f, x, p, lb, ub; alp_max=1f0, c=1f-4, ro=0.5f0)
    alp = alp_max
    g.=∇f(x)
    #f32max = floatmax(Float32)/10
   # g .= clamp.(g, -f32max, f32max)
    projected_gradient!(pg, x, g, lb, ub)
    
    # Compute directional derivative
    d = dot(pg, p)
    fx = f(x);
    # Find feasible alp
    while true
        x_new = clamp.(x .+ alp .* p, lb, ub)
        if f(x_new) <= fx + c * alp * d
            break
        end
        alp *= ro
        if alp < 1e-16
            break
            #error("Line search failed: step size too small")
        end
    end
    
    return alp
end


# Bounds: x1 ∈ [-0.5, 1.5], x2 ∈ [-1.0, 2.0]
x0 = [0.0, 0.0]
lb = [-0.5, -1.0]
ub = [0.5, 2.0]

x_opt = bounded_ncg(rosenbrock, ∇rosenbrock, x0, lb, ub, β_method="PR")
println("Optimal solution: x = ", x_opt)
println("Minimum value: f(x) = ", rosenbrock(x_opt))
println("Gradient at solution: ∇f(x) = ", ∇rosenbrock(x_opt))