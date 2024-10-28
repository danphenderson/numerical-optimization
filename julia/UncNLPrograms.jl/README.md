# UncNLPrograms
A set of high-dimensional, nonlinear, and unconstrained optimization problems implemented in
native Julia to test solvers using Automatic/Algorithmic Differentiation.

The UncNLPrograms is a subset of CUTEst, the latest evolution in the _Constrained and 
Unconstrained Testing Enviroment_. Each program defines an initial dimension, iterate, 
an objective function and its corresponding gradient.

The UncNLPrograms interface was inspired by the `NLPModels` package from
the JuliaSmoothOptimizers orginization.


## Interface
Function           | Description
-------------------|------------
`Programs`		   | Lists each selectable program and it's current dimension 
`SelectProgram`	   | Returns a reference to the specified program
`adjdim!`		   | Adjust the specified programs dimensions to the given order
`obj`			   | Evaluate the programs objective function at a point x
`grad`			   | Evaluate the programs gradient function at a point x
`objgrad`		   | Economical evaluation of the programs gradient and objective at a point x
`hessAD`		   | Returns the full Hessian matrix of an objective function at a point x (by Forward AD)

## Installation
Using Pkg.jl:
```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/danphenderson/UncNLPrograms.jl")
```  

Or from the REPL:
```julia
julia> ]
pkg> add "https://github.com/danphenderson/UncNLPrograms.jl"
``` 


## Simple UncNLPrograms Use Case
```julia
# Lists the name and dimension in the programs enviroment
Programs() 

# Store a reference to a program, as specified by the name given in Programs()
# programs are named as their origional .SIF file names
nlp = SelectProgram("SROSENBR")

# Returns the default itterate corresponding to the CUTEst problem in 20 dimensions
adjdim!(nlp, 20)

# Returns $f(x))$ at $x$ = nlp.$x0$, 
# where nlp.x0 is a reference to the default itterate
obj(nlp, nlp.x0) 

# Returns $\nabla x f(x))$ at $x$ = nlp.$x0$ 
grad(nlp, nlp.x0)

# Returns the tuple $(\nabla f(x),  \nabla f(x))$ at $x$ = nlp.$x0$ 
objgrad(nlp, nlp.x0)

# The jacobian of $\nabla f$ at $x$ = nlp.$x0$ using ForwardDiff.jl
hessAD(nlp, nlp.x0)
```

## Development Items
1. Determine a suitable home for UncNLPrograms in the Julia ecosystem and/or utilize an existing structure (i.e. build upon NLPModels, JuMP).
2. Optimize programs' implementation.
3. Check off each item in the Precompile warning list.
4. Overload Base.show(::UncProgram)
5. Confirm julia defined f and $\nabla^2 f(x)$ math the formula in Buckley's report.
