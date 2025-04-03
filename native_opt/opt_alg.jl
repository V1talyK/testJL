using LinearAlgebra

function mma(f, df, x0, bounds; max_iter=100, tol=1e-6, move=0.2f0)
    """
    Method of Moving Asymptotes (MMA) optimization algorithm in Julia.

    Parameters:
    f (function): Objective function to minimize.
    df (function): Gradient of the objective function.
    x0 (Vector): Initial guess for the design variables.
    bounds (Vector of Tuples): Lower and upper bounds for each design variable.
    max_iter (Int): Maximum number of iterations.
    tol (Float64): Convergence tolerance.
    move (Float64): Move limit for the asymptotes.

    Returns:
    x (Vector): Optimized design variables.
    f_val (Float64): Objective function value at the optimized point.
    """
    tp =eltype(x0)
    n = length(x0)
    x = copy(x0)
    L = [b[1] for b in bounds]  # Lower bounds
    U = [b[2] for b in bounds]  # Upper bounds
    alpha = zeros(tp, n)
    beta = zeros(tp, n)
    p = zeros(tp, n)
    q = zeros(tp, n)

    for k in 1:max_iter
        f_val = f(x)
        grad = df(x)
        println(k," ",f_val)
        println(k," ",grad[1:3])
        # Update asymptotes
        alpha = x - move * (U - L)
        beta = x + move * (U - L)

        # Update moving asymptotes
        for i in 1:n
            if grad[i] > 0
                p[i] = (U[i] - x[i])^2 * grad[i]
                q[i] = 0
            else
                p[i] = 0
                q[i] = -(x[i] - L[i])^2 * grad[i]
            end
        end

        # Update design variables
        x_new = zeros(tp, n)
        for i in 1:n
            tmp = 0f0;
            if !(iszero(p[i]) & iszero(q[i]))
                tmp = (p[i] - q[i]) / (2 * (p[i] + q[i])) * (beta[i] - alpha[i])
            end
            x_new[i] = (alpha[i] + beta[i]) / 2 - tmp
        end
        println(k," ",x_new[1:3], " " ,alpha[1:3]," ", beta[1:3]," ", p[1:3]," ", q[1:3], " ", L[1:3])
        # Apply bounds
        x_new = clamp.(x_new, L, U)

        # Check convergence
        if norm(x_new - x) < tol
            break
        end

        x = x_new
    end

    return x, f(x)
end

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    # Define objective function and its gradient
    f(x) = x[1]^2 + x[2]^2
    df(x) = [2 * x[1], 2 * x[2]]

    # Initial guess and bounds
    x0 = [5.0, 5.0]
    bounds = [(-10.0, 10.0), (-10.0, 10.0)]

    # Run MMA
    x_opt, f_opt = mma(f, df, x0, bounds)

    println("Optimized x: ", x_opt)
    println("Optimized f(x): ", f_opt)
end



function lrcB0n(XX::Matrix{Float32},
                YY::Array{Float32}, lb, ub; lf1=SEL, gf1=dSEL)
    #Линейная с ограничениями LBFGSB
    #minx = = similar(YY)
    eR = similar(YY)
    # z = similar(lb)
    # function Jf(a)
    #     mul!(YYc, XX, a)
    #     return sum(lf1.(YY, YYc))
    # end

    # function grad!(a)
    #     mul!(YYc, XX, a)
    #     #BLAS.gemv!('N', 1.0f0, XX, a, 0.0f0, YYc)
    #     eR .= gf1.(YY, YYc)
    #     z .= XX' * eR
    #     return z
    # end

    nx = size(XX, 2) # the dimension of the problem
    #bounds = collect(zip(lb, ub))
    # # opt.xtol_rel = 1e-5
    # # opt.maxeval = 100
    minx = ones(Float32, nx)
    grad = ones(Float32, nx)
    # x0 = clamp.(x0, lb, ub)

    #minx, minf = mma(Jf, grad!, x0, bounds; max_iter = 100)
    bvls!(minx, grad, eR, XX, YY, lb, ub)
    return XX * minx, minx
end



function interior_point(A, b, c; max_iter=100, tol=1e-6, mu=10.0, alpha=0.01)
    """
    Interior-Point Method for Linear Programming in Julia.

    Solves the problem:
        minimize c' * x
        subject to A * x = b, x >= 0

    Parameters:
    A (Matrix): Constraint matrix.
    b (Vector): Right-hand side vector.
    c (Vector): Cost vector.
    max_iter (Int): Maximum number of iterations.
    tol (Float64): Convergence tolerance.
    mu (Float64): Barrier parameter.
    alpha (Float64): Step size parameter.

    Returns:
    x (Vector): Optimized solution.
    """
    m, n = size(A)
    x = ones(n)  # Initial guess for x (must be strictly positive)
    λ = zeros(m)  # Lagrange multipliers for equality constraints
    s = ones(n)  # Slack variables (must be strictly positive)

    for iter in 1:max_iter
        # Residuals
        r_dual = A' * λ + s - c  # Dual feasibility
        r_pri = A * x - b        # Primal feasibility
        r_cent = x .* s          # Complementarity

        # Check convergence
        if norm(r_dual) < tol && norm(r_pri) < tol && norm(r_cent) < tol
            break
        end

        # Construct the KKT system
        X_inv = Diagonal(1.0 ./ x)
        S = Diagonal(s)
        K = [zeros(n, n) A' I;
             A zeros(m, m) zeros(m, n);
             S zeros(n, m) X_inv]

        # Right-hand side of the KKT system
        rhs = [-r_dual; -r_pri; -r_cent .+ mu]

        # Solve the KKT system
        Δ = K \ rhs
        Δx = Δ[1:n]
        Δλ = Δ[n+1:n+m]
        Δs = Δ[n+m+1:end]

        # Step size calculation
        α_primal = minimum([1.0; -alpha * x[Δx .< 0] ./ Δx[Δx .< 0]])
        α_dual = minimum([1.0; -alpha * s[Δs .< 0] ./ Δs[Δs .< 0]])

        # Update variables
        x += α_primal * Δx
        λ += α_dual * Δλ
        s += α_dual * Δs
    end

    return x
end

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    # Define problem data
    A = [1.0 1.0; 2.0 1.0; 3 0]  # Constraint matrix
    b = [4.0, 6.0, 1]          # Right-hand side
    c = [0.0, 0.0, 0.0]        # Cost vector

    # Solve using interior-point method
    x_opt = interior_point(A', b, c)

    println("Optimized x: ", x_opt)
end
