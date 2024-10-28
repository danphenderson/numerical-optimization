#    Problem : GROUP A
#    *********
#    This is a simplified version of problem NONDIA.
#
#    Original SIF Source:
#    G. Li,
#    "The secant/finite difference algorithm for solving sparse
#    nonlinear systems of equations",
#    SIAM Journal on Optimization, (to appear), 1990.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    LIARWHD.SIF classification  SUR2-AN-V-0
#
#    Number of variables (at least 2)
#
# Daniel Henderson, 08/2021

f = (x) -> begin
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        α = 2.0(x[i]^2 - x[1])
        γ = x[i]-1
        fx += α^2 + γ^2
    end
    return fx
end

g! = (g, x) -> begin
    for i in firstindex(x):lastindex(x)
        α = 2.0(x[i]^2 - x[1])
        γ = x[i]-1
        g[i] = 8.0x[i]α + 2.0γ
        g[1] -= 4.0α
    end
    return g
end

fg! = (g, x) -> begin
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        α = 2.0(x[i]^2 - x[1])
        γ = x[i]-1
        fx += α^2 + γ^2
        g[i] = 8.0x[i]α + 2.0γ
        g[1] -= 4.0α
    end
    return fx, g
end

init = (n::Int=5000) -> begin
	x0 = 4.0*ones(n)
    return n, x0
end

TestSet["LIARWHD"] = UncProgram("LIARWHD", f, g!, fg!, init)
