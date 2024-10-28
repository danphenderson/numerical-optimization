#    Problem : GROUP A
#    *********
#    The chained Rosenbrock function as given by Fletcher.
#
#    Original SIF Source:
#    The second problem given by
#    R. Fletcher,
#    "An optimal positive definite update for sparse Hessian matrices"
#    Numerical Analysis report NA/145, University of Dundee, 1992.
#
#    FLETCHR.SIF classification OUR2-AN-V-0
#
#    The Number of variables is variable.
#
# Daniel Henderson, 08/2021

f = (x) -> begin
    fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        γ = x[i+1] - x[i] + 1 - x[i]^2
        fx += 100*γ^2
    end
    return fx
end

g! = (g, x) -> begin
    for i in firstindex(x):lastindex(x)-1
        g[i] += 200*(1 + 2x[i])(-1+x[i]+x[i]^2-x[i+1])
        g[i+1] += 200*(1 - x[i] - x[i]^2 + x[i+1])
    end
    return g
end
 
fg! = (g, x) -> begin
    fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        γ = x[i+1] - x[i] + 1 - x[i]^2
        fx += 100*γ^2
        g[i] += 200*(-1 - 2x[i])*(1 - x[i] - x[i]^2 + x[1+i])
        g[i+1] += 200*(1 - x[i] - x[i]^2 + x[1+i])
    end
    return fx, g
end

init = (n::Int=1000) -> begin
	n < 2 && @warn("FLETCHCR: number of variables must be ≥ 2")
	n = max(n, 2)
    
    return n, zeros(n)
end

@warn "TODO: FLETCHCR, determine starting iterate - see whiteboard"
# TestSet["FLETCHCR"] = UncProgram("FLETCHCR", f, g!, fg!, init)