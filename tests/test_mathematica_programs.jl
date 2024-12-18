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

function read_input_data(file)
    """Reads input data test data of the form:
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

function test_mathematica_programs()
    """Tests the accuracy of our Mathematica test set of our variable
    unconstrained minimization problems.
    """
    # Read input data
    test_problems = read_input_data("public/fx0-dim10-mathematica.json")

    # Test each problem
    for test_problem in test_problems
        # Get the problem
        problem = eval(Meta.parse("ADNLPProblems.$test_problem[:name])(n=$n)"))
        
        expected_fx0 = obj(problem, problem.meta.x0)
        @test abs(expected_fx0 - test_problem.fx0)/expected_fx0 < 1e-12
        @test abs(expected_fx0 - test_problem.fx0) < 1e-8
    end
end
