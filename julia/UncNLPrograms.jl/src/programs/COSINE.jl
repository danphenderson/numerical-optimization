#    Problem : GROUP A
#    *********
#    Another function with nontrivial groups and repetitious elements.
#
#    Original SIF Source:
#    N. Gould, private communication.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    COSINE.SIF classification OUR2-AN-V-0
#
#    Number of variables (variable)

f = (x) -> begin
    fx = 0.0
    for i in 1:lastindex(x)-1
        α = -0.5x[i+1]+x[i]^2
        fx += cos(α)
    end
    return fx
end

g! = (g, x) -> begin
    for i in 1:lastindex(x)-1
        α = -0.5x[i+1]+x[i]^2
        g[i] -= 2.0sin(α)x[i]
        g[i+1] += 0.5sin(α)
    end
    return g
end
 
fg! = (g, x) -> begin
    fx = 0.0 
    for i in 1:lastindex(x)-1
        α = -0.5x[i+1]+x[i]^2
        fx += cos(α)
        g[i] -= 2.0sin(α)*x[i]
        g[i+1] += 0.5sin(α)
    end
    return fx, g
end

init = (n::Int=10000) -> begin
    n < 2 && @warn("COSINE: number of variables must be ≥ 2")
    return n, ones(n)
end

TestSet["COSINE"] = UncProgram("COSINE", f, g!, fg!, init)