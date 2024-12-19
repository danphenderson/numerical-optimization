using Plots
using LinearAlgebra
using Statistics
using OptimizationProblems, ADNLPModels, NLPModels
using DataFrames

# TODO: SampleBox should be an object

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

function build_sample_box(x0::Vector)
    return nothing
end


function pull_sample(box)
    """
    Pull a sample from the box
    """
    return x0 .+ (box .* randn(length(x0)))
end


function sample_around_x0(program, number_samples::Int=100)
    """
    Sample the objective function

    Note, we will want to collect samples around our minimizers within
    the step size of the minimizer.
    """
    box = build_sample_box(program.meta.x0)
    samples = [pull_sample(box) for _ in 1:number_samples]
    
    # Detrmine number of negative, positive, and singular eigenvalues
    # of the Hessian at the sampled point
    
    n_neg = 0
    n_pos = 0
    n_sing = 0

    for x in samples
        hess = hessian(program, x)
        eigvals = eigvals(hess)
        n_neg += sum(eigvals .< 0)
        n_pos += sum(eigvals .> 0)
        n_sing += sum(eigvals .== 0)
    end
    
    # Determine the percentage of negative, positive, and singular eigenvalues scaled by the magnitude of the Hessian?
    
    return n_neg, n_pos, n_sing
end

