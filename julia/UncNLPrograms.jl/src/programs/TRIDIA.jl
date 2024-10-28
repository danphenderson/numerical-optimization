#    Problem : GROUP A
#    *********
#	 Shanno's TRIDIA quadratic tridiagonal problem
#
#    Origonal SIF Source:problem 8 in
#    Ph.L. Toint,
#    "Test problems for partially separable optimization and results
#    for the routine PSPMIN",
#    Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
#
#    See also Buckley#40 (p.96)
#
#    TRIDIA.SIF classification QUR2-AN-V-0
#
#	 This problem is decomposed in n linear groups, the last n-1 of which
#    are 2 x 2 and singular.
#
#    Number of variables is variable

f = (x) -> begin
	fx   = (x[1]-1)^2
	for i in firstindex(x)+1:lastindex(x)
		item    = 2x[i]-x[i-1]
    	fx     += item^2*i
	end
    return fx
end

g! = (g, x) -> begin
	g[1] = 2(x[1]-1)
	for i in firstindex(x)+1:lastindex(x)
		item    = 2x[i]-x[i-1]
        g[i]   += 4item*i
	    g[i-1] -= 2item*i
	end
    return g
end

fg! = (g, x) -> begin
	fx   = (x[1]-1)^2
	g[1] = 2(x[1]-1)
	for i in firstindex(x)+1:lastindex(x)
		item    = 2x[i]-x[i-1]
    	fx     += item^2*i
        g[i]   += 4item*i
	    g[i-1] -= 2item*i
	end
    return fx, g
end

init = (n::Int=5000) -> begin
    return n, ones(n)
end

TestSet["TRIDIA"] = UncProgram("TRIDIA",f, g!, fg!, init)