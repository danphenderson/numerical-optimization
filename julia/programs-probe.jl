using Plots
using LinearAlgebra
using Optim, LineSearches
using OptimizationProblems, ADNLPModels, NLPModels
using DataFrames


"""
Module probes our test set by testing out various Conjugate Gradient methods

We are building on-top of ADNLPModels defined in OptimizationProblems.jl,
then we use Optim.jl to perform minimization on the models.
"""

function get_test_set(n::Int=40)
    """
    List of unconstrained scalable ADNLPModels within OptimizationProblems.jl
    """
    meta = OptimizationProblems.meta
    names_pb_vars = meta[
    (meta.variable_nvar .== true) .& (meta.ncon .== 0) .& (5 .<= meta.nvar .<= 100),
        [:nvar, :name]
    ]
    test_set_generator = (
        eval(Meta.parse("ADNLPProblems.$(pb[:name])(n=$n)")) for pb in eachrow(names_pb_vars)
    )
    return [p for p in test_set_generator]
end

function optim_options()
    """
    Using really stict conditions in low dimensions.

    Tacking earlier in routine may be obpuscated by extra
    iterations to obtain terminal convergence.
    """
    return Optim.Options(
        iterations = 100000,
        g_abstol = 1e-14,     
        f_abstol = 1e-14,     
        x_abstol = 1e-14,     
        store_trace = false, # Trace has a lot of useful stuff...
        show_trace = false,
        extended_trace = false,
    )
end

function bfgs_alg()
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

function qn_alg()
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

function run(problem::ADNLPModel, alg, opt)
    """
    Run the optimization algorithm on the problem
    """
    # Define objective and in-place gradient aligning with optim's interface.
    x0 = problem.meta.x0
    f(x) = obj(problem, x)
    g!(G, x) = grad!(problem, x, G) 

    # Run the optimization
    res = Optim.optimize(f, g!, x0, alg, opt)

    # Convert res object to a dataframe
    res = DataFrame(
        x = res.minimizer,
        fx = res.minimum,
        g_residual = norm(grad(problem, res.minimizer)),
        niter = res.iterations,
        neval_obj = res.f_calls,
        neval_grad = res.g_calls,
        trace = res.trace,
    )

    # Reset problem counters.
    NLPModels.reset!(problem)
    return res
end