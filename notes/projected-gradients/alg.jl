# Define a struct for bound constraints
struct Constraints
    lower::Vector{Float64}
    upper::Vector{Float64}
end

# (Optional) Define a struct for a quadratic objective to store Q and c, for clarity.
struct QuadraticObjective
    Q::Matrix{Float64}
    c::Vector{Float64}
end

# Make QuadraticObjective callable to compute q(x) = 0.5*x'*Q*x + c'*x.
function (obj::QuadraticObjective)(x::AbstractVector)
    return 0.5 * dot(x, obj.Q * x) + dot(obj.c, x)
end

# Helper to compute the gradient of q at x.
# If q is a QuadraticObjective (or has fields Q and c), use them directly; otherwise adapt as needed.
function _gradient(q, x::AbstractVector)
    if hasproperty(q, :Q) && hasproperty(q, :c)
        return q.Q * x + q.c
    else
        error("Cannot compute gradient: q must provide Hessian (Q) and linear term (c).")
    end
end

# Helper to compute Hessian times a vector (for quadratic q).
function _hess_mul(q, v::AbstractVector)
    if hasproperty(q, :Q)
        return q.Q * v
    else
        error("Cannot compute Hessian*vector: q must provide Hessian (Q).")
    end
end

# Check if a point x is feasible (satisfies l ≤ x ≤ u).
function _is_feasible(q, x::AbstractVector, constraints::Constraints)::Bool
    for i in eachindex(x)
        if x[i] < constraints.lower[i] || x[i] > constraints.upper[i]
            return false
        end
    end
    return true
end

# Check KKT conditions for the bound-constrained quadratic problem at point x.
function _is_kkt_satisfied(q, x::AbstractVector, constraints::Constraints)::Bool
    # Compute gradient at x
    local_g = _gradient(q, x)
    # Tolerance for near-zero gradient
    grad_tol = 1e-8
    for i in eachindex(x)
        local_l = constraints.lower[i]
        local_u = constraints.upper[i]
        if x[i] > local_l + 1e-12 && x[i] < local_u - 1e-12
            # x[i] is strictly inside bounds
            if abs(local_g[i]) > grad_tol
                return false  # interior point must have gradient ~ 0
            end
        elseif x[i] <= local_l + 1e-12
            # x[i] at (or very close to) lower bound
            if local_g[i] < -grad_tol
                return false  # gradient pointing significantly below feasible region
            end
        elseif x[i] >= local_u - 1e-12
            # x[i] at (or very close to) upper bound
            if local_g[i] > grad_tol
                return false  # gradient pointing significantly above feasible region
            end
        end
    end
    return true
end

# Compute the generalized Cauchy point for the quadratic function q at current point x_k under bounds.
function _find_cauchy_point(q, x_k::Vector{Float64}, constraints::Constraints)
    # Ensure starting point is feasible (project if not)
    local_x = similar(x_k)
    for i in eachindex(x_k)
        local_x[i] = min(max(x_k[i], constraints.lower[i]), constraints.upper[i])
    end
    # Initialize the set of fixed (active) constraints for this path
    n = length(local_x)
    fixed = falses(n)  # boolean array
    # Initial gradient and initial active set (constraints that are active with correct KKT sign)
    local_g = _gradient(q, local_x)
    for i in 1:n
        if (local_x[i] <= constraints.lower[i] + 1e-12 && local_g[i] >= 0) ||
           (local_x[i] >= constraints.upper[i] - 1e-12 && local_g[i] <= 0)
            fixed[i] = true
        end
    end

    # Follow the projected steepest descent path
    while true
        # Compute feasible descent direction d (projected negative gradient)
        local_d = -local_g
        for i in 1:n
            if fixed[i]
                # Do not move along directions that are fixed (active)
                local_d[i] = 0.0
            elseif local_x[i] <= constraints.lower[i] + 1e-12 && local_d[i] < 0
                # At lower bound and descent direction points outside feasible region
                local_d[i] = 0.0
                fixed[i] = true
            elseif local_x[i] >= constraints.upper[i] - 1e-12 && local_d[i] > 0
                # At upper bound and descent direction points outside feasible region
                local_d[i] = 0.0
                fixed[i] = true
            end
        end
        # If no feasible descent direction (all zero), then we are at a (projected) stationary point
        if all(local_d .== 0)
            break  # local_x remains as the Cauchy point (likely already KKT optimal)
        end

        # Determine step length to the nearest bound in direction d for each free component
        # and find the minimum positive step (first bound hit).
        t_min = Inf
        # We will collect indices that hit a bound at t_min (to handle ties)
        bound_hits = Int[]
        for i in 1:n
            if local_d[i] > 0
                # moving upward, check upper bound
                if constraints.upper[i] < Inf
                    local_t = (constraints.upper[i] - local_x[i]) / local_d[i]
                    if local_t < t_min - 1e-12
                        # Found a smaller step; reset hits list
                        t_min = local_t
                        bound_hits = [i]
                    elseif abs(local_t - t_min) < 1e-12
                        # Equal (within tolerance) to current min step, add to hits
                        push!(bound_hits, i)
                    end
                end
            elseif local_d[i] < 0
                # moving downward, check lower bound
                if constraints.lower[i] > -Inf
                    local_t = (constraints.lower[i] - local_x[i]) / local_d[i]  # local_d[i] < 0 gives positive t
                    if local_t < t_min - 1e-12
                        t_min = local_t
                        bound_hits = [i]
                    elseif abs(local_t - t_min) < 1e-12
                        push!(bound_hits, i)
                    end
                end
            end
        end

        # Compute directional derivative info for quadratic model
        # g^T d (linear term) and d^T Q d (quadratic term)
        local_dot_gd = dot(local_g, local_d)
        local_dQd = dot(local_d, _hess_mul(q, local_d))  # assuming q.Q is symmetric Hessian

        if local_dQd > 1e-12
            # Positive curvature: check if minimizer along d occurs before hitting a bound
            local_t_opt = - local_dot_gd / local_dQd
            if local_t_opt > 1e-12 && local_t_opt < t_min - 1e-12
                # A local minimum along the direction occurs before any bound is hit
                local_x .= local_x .+ local_t_opt * local_d
                # Ensure numerical feasibility (clip just in case of tiny overshoot)
                for i in 1:n
                    if local_x[i] < constraints.lower[i]
                        local_x[i] = constraints.lower[i]
                    elseif local_x[i] > constraints.upper[i]
                        local_x[i] = constraints.upper[i]
                    end
                end
                return local_x  # Found the generalized Cauchy point
            end
        end

        if t_min == Inf
            # No bound will be hit (direction unbounded within feasible region)
            # If curvature is non-positive, objective decreases without bound; we break to prevent infinite loop.
            # If curvature is positive, we'd have found t_opt above. So here either d^TQd <= 0 or domain is effectively unbounded.
            # We terminate by taking a "large" step (for safety, take 1.0 as step) and break.
            local_x .= local_x .+ (local_d .* 1.0)
            for i in 1:n
                if local_x[i] < constraints.lower[i]
                    local_x[i] = constraints.lower[i]
                elseif local_x[i] > constraints.upper[i]
                    local_x[i] = constraints.upper[i]
                end
            end
            break
        end

        # Take step t_min to the nearest bound and update x
        local_x .= local_x .+ t_min * local_d
        # Clamp those variables exactly to their bounds and mark as fixed
        for j in bound_hits
            if local_d[j] > 0
                local_x[j] = constraints.upper[j]
            elseif local_d[j] < 0
                local_x[j] = constraints.lower[j]
            end
            fixed[j] = true  # this bound is now active for the remainder of this path
        end

        # Compute new gradient after this step, then continue
        local_g = _gradient(q, local_x)
        # (The loop will recompute a new direction and possibly continue)
    end

    return local_x
end

# Solve the quadratic subproblem (16.74) given current point x (usually the Cauchy point),
# to find an improved feasible point x_plus with q(x_plus) <= q(x) (sufficient decrease).
function _solve_subproblem(q, x_current::Vector{Float64}, constraints::Constraints)
    # Start from the given feasible point (often the Cauchy point)
    x_plus = copy(x_current)
    n = length(x_plus)
    # Identify active and free sets at x_current
    active_idx = findall(i -> (x_current[i] <= constraints.lower[i] + 1e-12) || (x_current[i] >= constraints.upper[i] - 1e-12), 1:n)
    free_idx   = setdiff(1:n, active_idx)

    # If no free variables, x_current is already on all bounds (corner point)
    if isempty(free_idx)
        return x_plus  # nothing to optimize, return as is
    end

    # Pre-compute objective value at start for comparison
    q_val_start = q(x_current)
    # We will iteratively enforce bounds on free variables if needed
    max_iter = length(free_idx)  # at most free count iterations for adding bounds
    iter = 0
    while iter < max_iter
        iter += 1
        # Build sub-matrices for Hessian and vectors for gradient components
        # Q_ff (free-free block of Hessian), c_f (free part of linear term), and gradient at x_plus for free part
        # We can compute gradient and extract free part directly instead of building c_f:
        grad_full = _gradient(q, x_plus)
        grad_f = grad_full[free_idx]  # gradient on free variables at current x_plus

        # Form the Hessian block Q_ff and cross term Q_fa for free vs active
        # (if active set not empty; if empty, Q_fa term is nonexistent)
        Q_ff = q isa QuadraticObjective ? q.Q[free_idx, free_idx] : q.Q[free_idx, free_idx]  # assume QuadraticObjective for clarity
        # Compute right-hand side for optimum: -(Q_fa * x_active + c_f)
        # Instead of explicitly forming Q_fa, use gradient: grad_f = Q_ff*x_f + Q_fa*x_active + c_f,
        # so Q_ff*x_f + (Q_fa*x_active + c_f) = grad_f. Thus (Q_fa*x_active + c_f) = grad_f - Q_ff*x_f.
        # But since we want Q_ff*x_new_f = - (Q_fa*x_active + c_f),
        # we can compute b = Q_fa*x_active + c_f = Q_ff*x_plus_free - grad_f  (rearranging previous).
        # Actually: grad_f = Q_ff * x_free + Q_fa * x_active + c_f, so let b = Q_fa*x_active + c_f = grad_f - Q_ff * x_free.
        # Using current x_plus (which has some x_free, possibly updated by previous loop):
        if isempty(active_idx)
            # No active variables, then grad_f = Q_ff*x_free + c_f, so b = c_f = grad_f - Q_ff*x_free
            nothing
        end
        b = grad_f - Q_ff * x_plus[free_idx]  # this yields Q_fa*x_active + c_f

        # Solve for new free-variable values: Q_ff * x_free_new = -b
        # i.e., x_free_new = - Q_ff^{-1} * b
        # Use a linear solver (with care if Q_ff is indefinite or singular)
        x_free_new = similar(x_plus[free_idx])
        # Solve linear system. If Q_ff is symmetric indefinite, we use \ which will do an LU (or we could use LDL).
        # We'll attempt a solve and catch any errors due to singularity.
        solved = false
        try
            x_free_new .= -(Q_ff \ b)
            solved = true
        catch err
            # If solve fails (e.g., singular matrix), use gradient step as fallback.
            x_free_new .= x_plus[free_idx]  # start from current
            x_free_new .+= -grad_f  # one step of steepest descent in free subspace
            solved = true
        end

        # Update x_plus for free variables
        for (j_idx, var_idx) in enumerate(free_idx)
            x_plus[var_idx] = x_free_new[j_idx]
        end

        # Project any free variable that went out of bounds back to the bounds
        # and move it from free set to active set.
        new_active = Int[]
        for (j_idx, var_idx) in enumerate(free_idx)
            if x_plus[var_idx] < constraints.lower[var_idx] - 1e-12
                x_plus[var_idx] = constraints.lower[var_idx]
                push!(new_active, var_idx)
            elseif x_plus[var_idx] > constraints.upper[var_idx] + 1e-12
                x_plus[var_idx] = constraints.upper[var_idx]
                push!(new_active, var_idx)
            end
        end
        if !isempty(new_active)
            # Update active and free sets
            append!(active_idx, new_active)
            sort!(active_idx)  # keep sorted (not strictly necessary)
            free_idx = setdiff(1:n, active_idx)
            # If all variables became active, break
            if isempty(free_idx)
                break
            else
                # Continue loop to resolve with new active constraints
                continue
            end
        else
            # Found a candidate x_plus with all free vars within bounds
            break
        end
    end

    # Ensure feasibility (just a final clamp for safety)
    for i in 1:n
        x_plus[i] = min(max(x_plus[i], constraints.lower[i]), constraints.upper[i])
    end

    # Ensure sufficient decrease: accept x_plus only if q(x_plus) <= q(x_current)
    # (Use a small tolerance to account for numerical error)
    if q(x_plus) > q_val_start + 1e-12
        # If not improved (or worsened), fallback to the starting point (no progress)
        x_plus .= x_current
    end

    return x_plus
end

# Main projected gradient method for bound-constrained quadratic optimization (Algorithm 16.5).
function projected_gradients(q::Function, x0::Vector{Float64}, constraints::Constraints; max_iter::Int=1000)
    # Start from a feasible initial point (project if necessary)
    x = copy(x0)
    for i in eachindex(x)
        x[i] = min(max(x[i], constraints.lower[i]), constraints.upper[i])
    end

    # Main iteration loop
    for k in 1:max_iter
        # Check KKT optimality conditions
        if _is_kkt_satisfied(q, x, constraints)
            return x  # Found optimal (stationary) point
        end
        # Compute generalized Cauchy point
        x_c = _find_cauchy_point(q, x, constraints)
        # Perform subspace optimization (projected search step) from Cauchy point
        x_plus = _solve_subproblem(q, x_c, constraints)
        # Update iterate
        if all(x .== x_plus)
            # No change (could not improve), stop to avoid infinite loop
            return x_plus
        end
        x .= x_plus  # move to next iterate
    end

    # If max_iter reached without KKT satisfaction, return the last iterate
    return x
end
