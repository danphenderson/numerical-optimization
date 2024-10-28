#    Problem : GROUP A
#    *********
#	 A multi-dimensional variant of HUMPS, a two dimensional function
#    with a lot of humps. The density of humps increases with the
#    parameter ZETA, making the problem more difficult.
#	 
#	 The problem is nonconvex.
#
#    Origonal SIF Source:
#    Ph. Toint, private communication, 1997
#
#    GENHUMP.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable
#
# Daniel Henderson, 08/2021  
f = (x) -> begin
	fx = 0.0
    ζ = 20.0
    for i in firstindex(x):lastindex(x)-1
        fx += sin(ζ*x[i])^2*sin(ζ*x[i+1])^2+0.05*(x[i]^2+x[i+1]^2)
    end
    return fx
end

g! = (g, x) -> begin
    ζ = 20.0
    n = lastindex(x)

    g[1] = ζ*sin(2ζ*x[1])*sin(ζ*x[2])^2+0.1x[1]
    for i in 1:n-2
        g[i+1] = 0.2x[i+1]+40cos(20x[i+1])sin(20x[i])^2*sin(20x[i+1])+40cos(20x[i+1])sin(20x[i+1])sin(20x[i+2])^2
    end
    g[n] = 2ζ*cos(ζ*x[n])*sin(ζ*x[n-1])^2*sin(ζ*x[n])+0.1x[n]
    return g
end

fg! = (g, x) -> begin
    ζ = 20.0
    n = lastindex(x)
    fx = 0.0
    g[1] = ζ*sin(2ζ*x[1])*sin(ζ*x[2])^2+0.1x[1]
    for i in 1:n-2
        fx += sin(ζ*x[i])^2*sin(ζ*x[i+1])^2+0.05*(x[i]^2+x[i+1]^2)
        g[i+1] = 0.2x[i+1] + 40cos(20x[i+1])sin(20x[i])^2*sin(20x[i+1])+40cos(20x[i+1])sin(20x[i+1])sin(20x[i+2])^2
    end
    fx += sin(ζ*x[n-1])^2*sin(ζ*x[n])^2+0.05*(x[n-1]^2+x[n]^2)
    g[n] = 2ζ*cos(ζ*x[n])*sin(ζ*x[n-1])^2*sin(ζ*x[n])+0.1x[n]
    return g
end


init = (n::Int=5000) -> begin
	n < 2 && @warn("GENHUMPS: number of variables must be ≥ 2")
	n = max(n, 2)
    x0 = -506.2*ones(n)
    x0[1] = -506 
    return n, x0
end

TestSet["GENHUMPS"] = UncProgram("GENHUMPS", f, g!, fg!, init)