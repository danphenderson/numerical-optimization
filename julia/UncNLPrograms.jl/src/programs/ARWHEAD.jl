#    Problem : GROUP A
#    *********
#    A quartic problem whose Hessian is an arrow-head (downwards) with
#    diagonal central part and border-width of 1.
#
#    Origonal SIF Source:
#    A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
#    "Performance of a multifrontal scheme for partially separable
#    optimization",
#    Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    ARWHEAD.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable

f = (x) -> begin
    fx = 0.0
    xn_sqr = x[lastindex(x)]^2
    for i in 1:lastindex(x)-1
        fx += (x[i]^2 + xn_sqr)^2 + 3.0 - 4.0x[i]
    end
    return fx
end


g! = (g, x) -> begin
    gn = 0.0
    cxn = 4.0*x[lastindex(x)]
    xn_sqr = cxn^2/16.0
    for i in firstindex(x):lastindex(x)-1
        α = x[i]^2 + xn_sqr
        g[i] = 4.0x[i]*α - 4.0
        gn += cxn*α
    end
    g[lastindex(x)] = gn
    return g
end

fg! = (g, x) -> begin
    gn = fx = 0.0
    cxn = 4.0*x[lastindex(x)]
    xn_sqr = cxn^2/16.0
    for i in firstindex(x):lastindex(x)-1
        α = x[i]^2 + xn_sqr
        fx += (x[i]^2 + xn_sqr)^2 + 3.0 - 4.0*x[i]
        g[i] = 4.0*x[i]*α - 4.0
        gn += cxn*α
    end
    g[lastindex(x)] = gn
    return fx, g
end

init = (n::Int=5000) -> begin
	n < 2 && @warn("ARWHEAD: number of variables must be ≥ 2")
	n = max(n, 2)
    
    return n, ones(n)
end

TestSet["ARWHEAD"] = UncProgram("ARWHEAD", f, g!, fg!, init)
