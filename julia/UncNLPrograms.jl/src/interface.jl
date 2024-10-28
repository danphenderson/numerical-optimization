"""
	interface.jl

Defines the interface to each UncProgram in the UncNLPrograms enviroment.
The interface should feel familiar to interacting with NLPModels.jl,
apart of the JuliaSmoothOptimizers orginization.
"""

# macro from NLPModels.jl
macro lencheck(l, vars...)
    exprs = Expr[]
    for var in vars
      varname = string(var)
      push!(exprs, :(
        if length($(esc(var))) != $(esc(l))
          throw(DimensionError($varname, $(esc(l)), length($(esc(var)))))
        end
      ))
    end
    Expr(:block, exprs...)
end


"""
    obj

```math
f(x)
````
Evaluate the objective function at the point x.
"""
function obj(nlp::UncProgram, x::AbstractArray{<:Real})
    @lencheck nlp.n x
    return nlp.f(x)
end


"""
    grad

```math
∇f(x)
```
Evaluate the gradient of the objective function at a point x.
"""
function grad(nlp::UncProgram, x::AbstractArray{<:Real})
    @lencheck nlp.n x
    g = zeros(length(x))
    return nlp.g!(g, x)
end


"""
    objgrad

```math
f(x), ∇f(x)
```
Evaluate the gradient and it's objective function at a point x, 
by a iterating over the dimesions in an economical fashian.
"""
function objgrad(nlp::UncProgram, x::AbstractArray{<:Real})
    @lencheck nlp.n x
    g = zeros(length(x))
    return nlp.fg!(x, g)
end


"""
    hessAD

```math
∇²f(x)
```
Evaluate the full Hessian matrix at a point x by computing
the gradient of the Jacobian formula using ForwardDiff.jl,
a forward automatic differentation package. 
"""
function hessAD(nlp::UncProgram, x::AbstractArray{<:Real}) 
    @lencheck nlp.n x
    g = zeros(length(x))
    return ForwardDiff.jacobian(nlp.g!, g, x)
end


"""
    adjdim!

```math
f(x), ∇f(x)
```
Change the dimensions of a specified problem in the TestSet
"""
function adjdim!(nlp::UncProgram, n::Int)
    #@warn "adjdim!: This operation has not been tested to gaurentee similair output"
    n, x0 = nlp.init(n)
    nlp = UncProgram(nlp, n, x0)
    TestSet[nlp.name] = nlp
    return nlp
end


"""
    Programs

```math
f(x), ∇f(x)
```
A list of unconstrained nonlinear programming problems in the current testing enviroment. 
"""
function Programs()
    for nlp in values(TestSet)
        @printf "%s with dimension %d\n" nlp.name nlp.n
    end
    return keys(TestSet)
end # changes must be propogated to gHSTest


"""
    SelectProgram

```math
f(x), ∇f(x)
```
Returns an instance of UncProgram. 
"""
function SelectProgram(key::String; n=nothing)
    if key ∈ keys(TestSet)
        !isa(n, Nothing) && adjdim!(TestSet[key], n)
        return TestSet[key]
    end
    @warn "The program $key is not in the testing enviroment" 
end