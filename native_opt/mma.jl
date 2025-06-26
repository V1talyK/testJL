using LinearAlgebra

function mma(f, df, x0, lb, ub; max_iter=100, tol=1e-6, move=0.2)
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
    n = length(x0)
    x = copy(x0)
    alpha = zeros(n)
    beta = zeros(n)
    p = zeros(n)
    q = zeros(n)

    for k in 1:max_iter
        f_val = f(x)
        grad = df(x)

        # Update asymptotes
        alpha = x - move * (ub - lb)
        beta = x + move * (ub - lb)

        # Update moving asymptotes
        for i in 1:n
            if grad[i] > 0
                p[i] = (ub[i] - x[i])^2 * grad[i]
                q[i] = 0
            else
                p[i] = 0
                q[i] = -(x[i] - lb[i])^2 * grad[i]
            end
        end

        # Update design variables
        x_new = zeros(n)
        for i in 1:n
            x_new[i] = (alpha[i] + beta[i]) / 2 - (p[i] - q[i]) / (2 * (p[i] + q[i])) * (beta[i] - alpha[i])
        end

        # Apply bounds
        x_new = clamp.(x_new, lb, ub)

        # Check convergence
        if norm(x_new - x) < tol
            break
        end

        x .= x_new
        if any(isnan, x)
            println(k)
            wer
        end
    end

    return x, f(x)
end

# Example usage

    # Define objective function and its gradient
nx = 4
xx = rand(1000,nx)
af = [2.0, -4, 3.1, 5]
yy = xx*af

f(x) = sum(abs2, yy .-xx*x)
df(x) = -2.0*xx'*(yy .- xx*x)

    # Initial guess and bounds
x0 = ones(nx).+1
lb = -10*ones(nx)
ub = 10*ones(nx)

    # Run MMA
x_opt, f_opt = mma(f, df, x0, lb, ub; max_iter=1000, tol = 1e-8)

println("Optimized x: ", x_opt)
println("Optimized f(x): ", f_opt)
