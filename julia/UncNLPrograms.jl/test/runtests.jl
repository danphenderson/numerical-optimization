using UncNLPrograms
using Test, CUTEst, LinearAlgebra, Printf
import NLPModels as jo

@testset "Translation Testing" begin
    include("translation_tests.jl")
end

@testset "Interface Testing" begin
    include("interface_tests.jl")
end