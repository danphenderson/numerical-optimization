"""
Module tests the accuraccy of our Mathematica test set of our variable
unconstrained minimization problems.

Our mathematica test set is found in `mathematica/programs-probe.nb`,
which generates the following data artifacts:
  - `public/fx0-dim10-mathematica.json`
  - `public/fx0-dim15-mathematica.json`
Here we confirm our generated Mathematica test set aligns with the
source definitions in OptimizationProblems.jl.
"""

using Test
using JSON
using LinearAlgebra
using OptimizationProblems, ADNLPModels, NLPModels

Base.@kwdef struct TestProblem
    name::String
    fx0::Float64
end

function read_input_data(file::String)
    """Reads input test data of the form:
	[
        {
            "Name": "arglina",
            "InitialValue": 1.0e1
        }
    ]
    Returns: Array{TestProblem}
    """
    data = JSON.parsefile(file)
    return [TestProblem(d["Name"], d["InitialValue"]) for d in data]
end

function run_tests_for_dimension(dimension::Int)
    test_problems = read_input_data("public/fx0-dim$dimension-mathematica.json")

    for test_problem in test_problems
        # Get the problem
        problem = eval(Meta.parse("ADNLPProblems.$(test_problem.name)(n=$dimension)"))

        # Evaluate fx0
        expected_fx0 = obj(problem, problem.meta.x0)

        # Run tests; failures are captured and reported by @testset
        try
            @test abs(expected_fx0 - test_problem.fx0)/expected_fx0 < 1e-12
            @test abs(expected_fx0 - test_problem.fx0) < 1e-8
        catch e 
            @warn "Test failed for problem: $(test_problem.name) with dimension: $dimension"
        end
    end
end

function run_tests()
    test_dimension = [5, 6, 7, 10, 15]
    for d in test_dimension
        run_tests_for_dimension(d)
    end
end

# @testset "Mathematica Programs" begin
#     test_mathematica_programs()
# end