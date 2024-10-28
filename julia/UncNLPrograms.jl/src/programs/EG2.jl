#    Problem : GROUP A
#    *********
#    TA simple nonlinear problem given as an example in Section 1.2.4 of
#    the LANCELOT Manual.
#    The problem is non convex and has several local minima.
#
#    Original SIF Source:
#    A.R. Conn, N. Gould and Ph.L. Toint,
#    "LANCELOT, A Fortran Package for Large-Scale Nonlinear Optimization
#    (Release A)"
#    Springer Verlag, 1992.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    EG2.SIF classification OUR2-AN-1000-0
#
#    Number of variables (at least 2)
#
# Daniel Henderson, 08/2021

f = (x) -> begin
    fx = 0.0 
    α = x[1] - 1.0
    for i in firstindex(x):lastindex(x)-1
        γ = α + x[i]^2
        fx += sin(γ)
    end
    return fx
end


g! = (g, x) -> begin
    α = x[1] - 1
    xn = x[lastindex(x)] 
    g[lastindex(x)] = cos(xn^2)xn
    for i in firstindex(x):lastindex(x)-1
        γ = α + x[i]^2
        g[1] += cos(γ)
        g[i] += 2.0cos(γ)x[i]
    end
    return g
end

fg! = (g, x) -> begin
    fx = 0.0
    α = x[1] - 1.0
    xn = x[lastindex(x)] 
    g[lastindex(x)] = cos(xn^2)xn
    for i in firstindex(x):lastindex(x)-1
        γ = α + x[i]^2
        fx += sin(γ)
        g[1] += cos(γ)
        g[i] += 2.0cos(γ)x[i]
    end
    return fx, g
end

init = (n::Int=1000) -> begin
    n < 2 && @warn("EG2: number of variables must be ≥ 2")
    n = max(n, 2)
    
    return n, zeros(n)
end 

#TestSet["EG2"] = UncProgram("EG2", f, g!, fg!, init)