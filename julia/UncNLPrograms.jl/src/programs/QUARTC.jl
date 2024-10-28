#    Problem : GROUP A
#    *********	 
#    A simple quartic function.
#
#    Origonal SIF Source: problem 157 (p. 87) in
#    A.R. Buckley,
#   "Test functions for unconstrained minimization",
#    TR 1989CS-3, Mathematics, statistics and computing centre,
#    Dalhousie University, Halifax (CDN), 1989.
#
#    QUARTC.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable
#
# Daniel Henderson, 08/2021   

f = x -> begin
    return sum((x[i] - i)^4 for i in 1:lastindex(x))	
end

g! = (g, x) -> begin
	for i in 1:lastindex(x)
		g[i] = 4(x[i] - i)^3
	end
    return g
end

fg! = (g, x) -> begin
	fx = 0.0
	for i in 1:lastindex(x)
		fx  += (x[i] - i)^4
		g[i] = 4(x[i] - i)^3
	end
    return fx, g
end

init = (n::Int=5000) -> begin
    return n, 2.0ones(n)
end

TestSet["QUARTC"] = UncProgram("QUARTC",  f, g!, fg!, init)