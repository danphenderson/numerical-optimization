using Arpack, Optim, LineSearches, ThreadsX, LDLFactorizations
using StaticArrays, ForwardDiff, OrdinaryDiffEq

include("common.jl")

function get_optim_options(; gtol = 1e-8, maxiter = 50_000, trace = false)
    return Optim.Options(
        iterations   = maxiter,   # hard stop
        g_abstol     = gtol,      # practical accuracy
        store_trace  = trace,
        show_trace   = trace,
        extended_trace = trace,
        show_warnings = true,
    )
end

function newton_trust_region()
    """
    Defines the Newton Trust Region algorithm for the optimization.
    
    TODO: Add citation to Nocedal and Wright why we use
    the following parameter values.
    """
    return NewtonTrustRegion(; 
        initial_delta = 1.0,
        delta_hat = 100.0,
        eta = 0.1,
        rho_lower = 0.25,
        rho_upper = 0.75
)
end

function newton_linesearch()
    """
    Defines the Newton Line Search algorithm for the optimization.
    """
    return Newton(; alphaguess = LineSearches.InitialStatic(),
    linesearch = LineSearches.HagerZhang())
end

function newton_static()
    """
    Defines the Newton algorithm for the optimization.
    This is the most basic version of the algorithm.
    It uses a static line search and a static step size,
    which matches theorems used where knowing convexity
    allows us to asset iterations bound.
    """
    return Newton(linesearch  = LineSearches.Static(),  # α ≡ 1
        alphaguess = LineSearches.InitialStatic())
end

function gradient_descent_linesearch()
    """
    Defines the Quasi-Newton algorithm for the optimization.

    P is our H. Currently P = \nabla^2 f(x) = I and we fallback
    to gradient descent. 
    
    TODO: Accept P as an argument.
    """
    GradientDescent(; 
        alphaguess = LineSearches.InitialHagerZhang(),
        linesearch = LineSearches.HagerZhang(),
        P = nothing,
        precondprep = (P, x) -> nothing
    )
end



function eigs_hess(problem::ADNLPModel, x_cp::Vector)
    H = hess(problem, x_cp)
    λ, _ = eigs(H, nev=problem.meta.nvar - 1, maxiter=10000, which=:LM)
    λmin, _ = eigs(H, nev=1, maxiter=10000, which=:SR)
    push!(λ, pop!(λmin))
    return λ
end


function maximize(problem::ADNLPModel, sample::Vector)
    """
    Run the optimization algorithm on the problem
    """
    # Define objective and in-place gradient aligning with optim's interface.
    x0 = sample
    f(x) = -obj(problem, x)
    g!(G, x) = -1 .* grad!(problem, x, G)
    h!(H, x) = copyto!(H, hess(problem, x))

    # Run the optimization
    res = Optim.optimize(
        f, 
        g!, 
        h!,
        x0, 
        newton_trust_region(),
        get_optim_options()
    )

    # Reset problem counters.
    NLPModels.reset!(problem)
    return res
end


function minimize(problem::ADNLPModel, sample::Vector)
    """
    Run the optimization algorithm on the problem
    """
    # Define objective and in-place gradient aligning with optim's interface.
    x0 = sample
    f(x) = obj(problem, x)
    g!(G, x) = grad!(problem, x, G)
    h!(H, x) = copyto!(H, hess(problem, x))

    # Run the optimization
    res = Optim.optimize(
        f, 
        g!, 
        h!,
        x0, 
        newton_trust_region(),
        get_optim_options()
    )

    # Reset problem counters.
    NLPModels.reset!(problem)
    return res
end


function simulate(problem::ADNLPModel)
    box = build_sample_box(problem)
    df = DataFrame(
        "initial_objective" => Vector{Float64}(),
        "final_objective" => Vector{Float64}(), 
        "g_residual" => Vector{Float64}(),
        "is_saddle" => Vector{Bool}(),
        "pos_curvature_directions" => Vector{Int}(),
        "neg_curvature_directions" => Vector{Int}(),
        "zero_curvature_directions" => Vector{Int}(),
        "max_lambda" => Vector{Float64}(),
        "min_lambda" => Vector{Float64}(),
        "median_lambda" => Vector{Float64}(),
        "iterations" => Vector{Int}(),
        "critical_point" => Vector{Vector{Float64}}(),  
    )
    for _ in 1:100
        sample = pull_sample(problem, box)
        fx0 = obj(problem, sample)
        try
            res = run_sample(problem, sample)
            λ = eigs_hess(problem, res.minimizer)
            push!(df, (
                initial_objective = fx0,
                final_objective = res.minimum,
                g_residual = res.g_residual,
                is_saddle = λ[end] < 0,
                pos_curvature_directions = sum(λ .> 0),
                neg_curvature_directions = sum(λ .< 0),
                zero_curvature_directions = sum(abs.(λ) .< 1e-12),
                max_lambda = maximum(λ),
                min_lambda = minimum(λ),
                median_lambda = median(λ),
                iterations = res.iterations,
                critical_point = res.minimizer,
            ))
        catch
            continue
        end
    end
    return df
end

function unique_critical_points(df::DataFrame)
    """
    Compares the critical points locations to determine
    the number of unique critical points.
    """
    critical_points = df.critical_point
    critical_points = [round.(x, digits=4, base=2) for x in critical_points] # HACK d
    return length(unique(critical_points))
end

function run_all()
    test_set = get_test_set()
    mkpath("public/saddles/dim-40")
    for pb in test_set
        problem = get_problem(pb, 40)
        df = run(problem)
        CSV.write("public/saddles/dim-40/$(pb).csv", df)
        println("$(pb) has $(unique_critical_points(df)) unique critical points.")
        println("   $(pb) total saddles $(sum(df.is_saddle))")
    end
end
