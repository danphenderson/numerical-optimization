using UncNLPrograms
using Test, CUTEst, LinearAlgebra, Printf
import NLPModels as jo

# Tests translation of CUTE problems against CUTEst
for nlp ∈ values(UncNLPrograms.TestSet)
    @printf "\nTesting %s to CUTEst.jl Model at default iterate:\n" nlp.name
    qt_nlp = CUTEstModel(nlp.name) 
    qt_x0  = qt_nlp.meta.x0
    x0 = nlp.x0
    
    # set test points in a random proximity of the default itterate, both +/=
    t1 = x0 - rand(nlp.n)
    t2 = x0 - rand(nlp.n)
    t3 = x0 - rand(nlp.n)

    t4 = x0 + rand(nlp.n)
    t5 = x0 + rand(nlp.n)
    t6 = x0 + rand(nlp.n)
    try
        # test similarity of default iterates
        @test obj(nlp, x0) ≈ jo.obj(qt_nlp, qt_x0)
        @test grad(nlp, x0) ≈ jo.grad(qt_nlp, qt_x0)
        @test hessAD(nlp, x0) ≈ jo.hess(qt_nlp, qt_x0)
        
        # test at 6 randomly generated points
        @test obj(nlp, t1) ≈ jo.obj(qt_nlp, t1)
        @test grad(nlp, t1) ≈ jo.grad(qt_nlp, t1)
        @test hessAD(nlp, t1) ≈ jo.hess(qt_nlp, t1)

        @test obj(nlp, t2) ≈ jo.obj(qt_nlp, t2)
        @test grad(nlp, t2) ≈ jo.grad(qt_nlp, t2)
        @test hessAD(nlp, t2) ≈ jo.hess(qt_nlp, t2)
        @test obj(nlp, t3) ≈ jo.obj(qt_nlp, t3)
        @test grad(nlp, t3) ≈ jo.grad(qt_nlp, t3)
        @test hessAD(nlp, t3) ≈ jo.hess(qt_nlp, t3)

        @test obj(nlp, t4) ≈ jo.obj(qt_nlp, t4)
        @test grad(nlp, t4) ≈ jo.grad(qt_nlp, t4)
        @test hessAD(nlp, t4) ≈ jo.hess(qt_nlp, t4)

        @test obj(nlp, t5) ≈ jo.obj(qt_nlp, t5)
        @test grad(nlp, t5) ≈ jo.grad(qt_nlp, t5)
        @test hessAD(nlp, t5) ≈ jo.hess(qt_nlp, t5)

        @test obj(nlp, t6) ≈ jo.obj(qt_nlp, t6)
        @test grad(nlp, t6) ≈ jo.grad(qt_nlp, t6)
        @test hessAD(nlp, t6) ≈ jo.hess(qt_nlp, t6)
    catch e
        println(e)
    end

    finalize(qt_nlp)
end