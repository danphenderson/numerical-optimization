#    Problem :
#    *********
#	 The separable extension of Rosenbrock's function.
#
#    Origonal SIF Source:problem 21 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#
#    SROSENBR.SIF classification SUR2-AN-V-0
#
#    Number of variables is variable 
#
# Daniel Henderson, 08/2021


f = (x) -> begin
	fx = 0.0
	for i in firstindex(x):2:lastindex(x)
		t1     = 1 - x[i]
		t2     = 10(x[i+1] - x[i]^2)
		fx    += t1^2 + t2^2
	end
    return fx
end

g! = (g, x) -> begin
	for i in firstindex(x):2:lastindex(x)
		t1     = 1 - x[i]
		t2     = 10(x[i+1] - x[i]^2)
		g[i+1] = 20t2
		g[i]   = -2(x[i] * g[i+1] + t1)
	end
    return g
end

fg! = (g, x) -> begin
	fx = 0.0
	for i in firstindex(x):2:lastindex(x)
		t1     = 1 - x[i]
		t2     = 10(x[i+1] - x[i]^2)
		g[i+1] = 20t2
		g[i]   = -2(x[i] * g[i+1] + t1)
		fx    += t1^2 + t2^2
	end
    return fx, g
end

init = (n::Int=5000) -> begin
	mod(n, 2) > 0 && @warn "SROSENBR: number of variables must be even" 
	q = max(1, div(n, 2))
	n = 2q

    x0 = [j % 2 == 1 ? 1.2 : 1.0 for j in 1:n]
    return n, x0
end

TestSet["SROSENBR"] = UncProgram("SROSENBR",  f, g!, fg!, init)