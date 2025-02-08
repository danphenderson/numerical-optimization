using CSV, DataFrames, OptimizationProblems, ADNLPModels, NLPModels
using Random, Arpack, Optim, LineSearches
using LinearAlgebra, Statistics, Distributions

Random.seed!(1234)


function get_test_set(n::Int=40)
    """
    List of unconstrained scalable ADNLPModels within OptimizationProblems.jl
    """
    # meta = OptimizationProblems.meta
    # names_pb_vars = meta[
    # (meta.variable_nvar .== true) .& (meta.ncon .== 0) .& (5 .<= meta.nvar .<= 100),
    #     [:nvar, :name]
    # ]
    # test_set_generator = (
    #     eval(Meta.parse("ADNLPProblems.$(pb[:name])(n=$n)")) for pb in eachrow(names_pb_vars)
    # )
    # return [p for p in test_set_generator]
    return [
	"genrose",
	"arglina",
	"freuroth",
	"eg2",
	"cosine",
	"arglinb",
	"arglinc",
	"argtrig",
	"arwhead",
	"bdqrtic",
	"brownal",
	"broyden3d",
	"chnrosnb_mod",
	"cragglvy",
	"cragglvy2",
	"curly10",
	"curly10",
	"curly20",
	"curly30",
	"dixon3dq",
	"dqdrtic",
	"dqrtic",
	"edensch",
	"engval1",
	"errinros_mod",
	"extrosnb",
	"fletcbv2",
	"fletcbv3_mod",
	"fletchcr",
	"genhumps",
	"genrose_nash",
	"indef_mod",
	"integreq",
	"liarwhd",
	"morebv",
	"noncvxu2",
	"noncvxun",
	"nondia",
	"nondquar",
	"penalty1",
	"penalty2",
	"penalty3",
	"power",
	"quartc",
	"sbrybnd",
	"schmvett",
	"scosine",
	"sinquad",
	"tointgss",
	"tquartic",
	"tridia",
	"vardim"];
end

function get_problem(name::String, n::Int=40)
    """
    Get a problem by name
    """
    return eval(Meta.parse("ADNLPProblems.$(name)(n=$n)"))
end

function get_optim_options()
    """
    Using really stict conditions in low dimensions.

    Tacking earlier in routine may be obpuscated by extra
    iterations to obtain terminal convergence.
    """
    return Optim.Options(
        iterations = 10000000,
        g_abstol = eps(),       
        store_trace = false, # Trace has a lot of useful stuff...
        show_trace = false,
        extended_trace = false,
    )
end

function build_sample_box(problem::ADNLPModel)
    """ 
    Box centered around the minimizer
    """
    x0 = problem.meta.x0
    scale_vector = 2 .* abs.(x0)
    scale_vector[scale_vector .<= 1.0] .= 1.0
    return scale_vector
end

function pull_sample(problem, box::Vector)
    """
    Pulls a sample uniformly from the box box surrounding x0
    """
    return rand.(Uniform.(-box, box))
end

function bfgs_linesearch()
    """
    Define the algorithm for the optimization
    """
    return BFGS(;
        alphaguess = LineSearches.InitialStatic(),
        linesearch = LineSearches.HagerZhang(),
        initial_invH = x -> Matrix{eltype(x)}(I, length(x), length(x)),
        manifold = Flat(),
    )
end

function gradiant_descent_linesearch()
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

function newton_trust_region()
    """
    Defines the Newton Trust Region algorithm for the optimization.
    """
    return NewtonTrustRegion(; initial_delta = 1.0,
        delta_hat = 100.0,
        eta = 0.1,
        rho_lower = 0.25,
        rho_upper = 0.75)
end

function eigs_hess(problem::ADNLPModel, x_cp::Vector)
    H = hess(problem, x_cp)
    λ, _ = eigs(H, nev=problem.meta.nvar - 1, maxiter=10000, which=:LM)
    λmin, _ = eigs(H, nev=1, maxiter=10000, which=:SR)
    push!(λ, pop!(λmin))
    return λ
end

function run_sample(problem::ADNLPModel, sample::Vector)
    """
    Run the optimization algorithm on the problem
    """
    # Define objective and in-place gradient aligning with optim's interface.
    x0 = sample
    f(x) = obj(problem, x)
    g!(G, x) = grad!(problem, x, G)
    h!(H, x) = hess!(problem, x, H)

    # Run the optimization
    res = Optim.optimize(
        f, 
        g!, 
        x0, 
        bfgs_linesearch(),
        get_optim_options()
    )

    # Reset problem counters.
    NLPModels.reset!(problem)
    return res
end

function run(problem)
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
    for _ in 1:1000
        sample = pull_sample(problem, box)
        fx0 = obj(problem, sample)
        try
            res = run_sample(problem, sample)
            if !res.f_converged
                continue
            end
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
