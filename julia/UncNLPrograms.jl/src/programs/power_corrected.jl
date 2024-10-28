#    Problem : GROUP C
#    *********
#    This is a Power problem by Oren, the corrected SIF version
# 	 of the POWER problem to reflect Oren's formula.
#
#    Origonal SIF Source:
#    S.S. Oren,
#    Self-scaling variable metric algorithms,
#    Part II: implementation and experiments"
#    Management Science 20(5):863-874, 1974.
#
#    POWER.SIF classification OUR2-AN-V-0
#
# Daniel Henderson, 08/2021
@warn "TODO: verify power_corrected.jl formula matches Oren's specification"  

f = (x) -> begin
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        fx += (x[i]i)^2
    end
    return fx
end

g! = (g, x) -> begin
    for i in firstindex(x):lastindex(x)
        g[i] = 2x[i]i^2
    end
    return g
end

fg! = (g, x) -> begin
    α = 0.0
    for i in firstindex(x):lastindex(x)
        α += x[i]^2*i
    end
    for i in firstindex(x):lastindex(x)
        g[i] = 4α*i*x[i]
    end
    return fx, g
end

init = (n::Int=10000) -> begin
    return n, ones(n)
end

@warn "power_corrected: Currently excluded, but maybe include after testing?"
#TestSet["power_corrected"] = UncProgram("power_corrected", f, g!, fg!, init)