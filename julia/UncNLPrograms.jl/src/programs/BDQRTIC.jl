#    Problem : GROUP A
#    *********
#    This problem is quartic and has a banded Hessian with bandwidth = 9.
#
#    Origonal SIF Source:
#    Problem 61 in
#    A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
#    "Performance of a multifrontal scheme for partially separable
#    optimization",
#    Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
#
#    Implementation translated from Source:
#    http://eprints.tsu.ge/234/14/Tests%20collection-K-F.pdf
#
#    BDQRTIC.SIF classification SUR2-AN-V-0
#
#    Number of variables (variable)

f = (x) -> begin
    fx = 0.0
    xn_sqr = x[lastindex(x)]^2
    for i in firstindex(x):lastindex(x)-4
        α = 3.0 - 4.0x[i]
        γ = x[i]^2 + 2.0x[i+1]^2 + 3.0x[i+2]^2 + 4.0x[i+3]^2 + 5.0xn_sqr
        fx += α*γ + γ^2
    end
    return fx
end

g! = (g, x) -> begin
    xn = x[lastindex(x)]
    cxn_sqr = 5.0xn^2
    gn = 0.0
    for i in firstindex(x):lastindex(x)-4
        α = 3.0 - 4.0x[i]
        γ = x[i]^2 + 2.0x[i+1]^2 + 3.0x[i+2]^2 + 4.0x[i+3]^2 + cxn_sqr
        g[i] += -8.0α + 4.0γ*x[i];
        g[i+1] += 8.0γ*x[i+1];
        g[i+2] += 12.0γ*x[i+2];
        g[i+3] += 16.0γ*x[i+3];
        gn += 20.0γ*xn;
    end
    g[lastindex(x)] = gn
    return g
end

fg! = (g, x) -> begin
    fx = 0.0 
    xn = x[lastindex(x)]
    cxn_sqr = xn^2
    gn = 0.0
    for i in firstindex(x):lastindex(x)-4
        α = 3.0 - 4.0x[i]
        γ = x[i]^2 + 2.0x[i+1]^2 + 3.0x[i+2]^2 + 4.0x[i+3]^2 + cxn_sqr
        fx += α*γ + γ^2
        g[i]   += -8.0α + 4.0γ*x[i];
        g[i+1] += 8.0γ*x[i+1];
        g[i+2] += 12.0γ*x[i+2];
        g[i+3] += 16.0γ*x[i+3];
        gn += 20.0γ*xn;
    end
    g[lastindex(x)] = gn
    return fx, g
end

init = (n::Int=5000) -> begin
	n < 5 && @warn("BDQRTIC: number of variables must be ≥ 5")
	n = max(n, 5)

    return n, ones(n)
end

@warn "TODO: Debug BDQRTIC"
#TestSet["BDQRTIC"] = UncProgram("BDQRTIC", f, g!, fg!, init)
