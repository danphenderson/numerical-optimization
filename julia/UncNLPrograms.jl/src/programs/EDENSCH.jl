#    Problem : GROUP A
#    *********
#    The extended Dennis and Schnabel problem, as defined by Li.
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
#    EDENSCH.SIF classification  OUR2-AN-V-0
#
#    Number of variables (at least 2)
#
# Daniel Henderson, 08/2021

f = (x) -> begin
    fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        item1 = x[i] - 2
        item2 = x[i]x[i+1] - 2x[i+1]
        item3 = x[i+1] + 1
        fx += 16 + item1^4 + item2^2 + item3^2
    end
    return fx
end

g! = (g, x) -> begin
    for i in firstindex(x):lastindex(x)-1
        item1 = x[i] - 2
        item2 = x[i]*x[i+1] - 2x[i+1]
        item3 = x[i+1] + 1
        g[i] += 4item1^3 + 2x[i+1]item2
        g[i+1] += 2item2*(x[i] - 2) + 2item3
    end
    return g
end

fg! = (g, x) -> begin
    fx = 0.0
    for i in firstindex(x):lastindex(x)-1
        item1 = x[i] - 2
        item2 = x[i]*x[i+1] - 2x[i+1]
        item3 = x[i+1] + 1
        fx += 16 + item1^4 + item2^2 + item3^2
        g[i] += 4item1^3 + 2x[i+1]item2
        g[i+1] += 2item2*(x[i] - 2) + 2item3
    end
    return fx, g
end

init = (n::Int=2000) -> begin

    return n, 8.0*ones(n)
end

@warn "TODO: Debug EDENSCH"
# TestSet["EDENSCH"] = UncProgram("EDENSCH", f, g!, fg!, init)