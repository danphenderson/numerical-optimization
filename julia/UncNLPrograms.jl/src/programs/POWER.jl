#    Problem : GROUP C
#    *********
#    This is an error of a Power problem by Oren.
#
#    Origonal SIF Source:
#    S.S. Oren,
#    Self-scaling variable metric algorithms,
#    Part II: implementation and experiments"
#    Management Science 20(5):863-874, 1974.
#
#    The CUTE implementation specification of the source 
#    formula is not correct. 
#
#    POWER.SIF classification OUR2-AN-V-0
#
# Daniel Henderson, 08/2021  

f = (x) -> begin
    fx = 0.0
    for i in firstindex(x):lastindex(x)
        fx += x[i]^2*i
    end
    return fx^2
end

g! = (g, x) -> begin
    α = 0.0
    for i in firstindex(x):lastindex(x)
        α += x[i]^2*i
    end
    for i in firstindex(x):lastindex(x)
        g[i] = 4α*i*x[i]
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
    return α^2, g
end

init = (n::Int=10000) -> begin
    return n, ones(n)
end

TestSet["POWER"] = UncProgram("POWER", f, g!, fg!, init)