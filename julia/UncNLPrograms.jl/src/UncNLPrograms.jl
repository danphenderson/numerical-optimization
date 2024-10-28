"""
	UncNLPrograms

A subset of high-dimensional nonlinear C² CUTEr/st problems in native Julia.
Each UncProgram in the UncNLPrograms specifies an analytical representation of a nonlinear
program f and the corresponding gradient ∇f. Additionally, each program has a method for computing the
(f, ∇f) over an economical iteration of the dimensions.

Future development of the UncNLPrograms should focus on the specification of a standardized
native julia constrained and unconstrained optimization enviroment. UncNLPrograms.jl contains
a small translation of CUTE, which is the standard suite for performing numerical 
expeirments in optimization research. The translation exists because a new Standard is needed,
as for our case to test accelerated schemes utilizing the SIMD parallel nature of Automatic 
Differentation. The Julia Langague offers the needed flexibilty in Automatic Differentation
computations, primarly through it's multiple dispatch design which enables operator overloading.
Furthermore, Julia has a very eleqant treatment of type,  allowing for a  and

There are many attempts at creating an underlying data structure to hold Standard Form
Programs. The need for a common underlying data structure is to facilitate the definition
of Optimization solving routines. The most supported modeling structure is JuMP, who
acknowledges the diversity and need for a standard program model in the ecosystem.

The UncNLPrograms remains agnostic to existing Standard Form programming problem
models, as none support flexible automatic differentation. Rather, it serves as the
testing set for an AD based quasi-newton scheme. See QuasiNewtonHS.jl
"""
module UncNLPrograms

using LinearAlgebra, Printf, ForwardDiff


"""
    UncProgram

A base parent type of each unconstrained non-linear program
"""
struct UncProgram
    name::String
    f::Function
    g!::Function
    fg!::Function
    init::Function
    n::Int
    x0::AbstractArray
    function UncProgram(name, f, g!, fg!, init)
        n, x0 = init()
        new(name, f, g!, fg!, init, n, x0)
    end 
    function UncProgram(nlp::UncProgram, n, x0)
        new(nlp.name, nlp.f, nlp.g!, nlp.fg!, nlp.init, n, x0)
    end 
end

include("interface.jl")

"""
    UncNLPrograms.TestSet

A dictionary mapping the problem name it's corresponding UncProgram.
"""
TestSet = Dict{String, UncProgram}()
for p in readdir(joinpath(@__DIR__, "programs"))
    include(joinpath("programs", p))
end



export obj, grad, objgrad, adjdim!, hessAD, Programs, SelectProgram

end
