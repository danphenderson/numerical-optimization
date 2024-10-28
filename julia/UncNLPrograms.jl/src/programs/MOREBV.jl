#    Problem :
#    *********	 
#
#    The Boundary Value problem.
#    This is the nonlinear least-squares version without fixed variables.
#  
#
#    Source:  problem 28 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#	
#    Implementation translated from:
#    See also Buckley#17 (p. 94):
#	
#	
#    SIF input: Ph. Toint, Dec 1989 and Nick Gould, Oct 1992.
#
#    classification SUR2-MN-V-0
#    The number of variables is N. 
#
# Daniel Henderson, 08/2021   

f = x -> begin
	n = lastindex(x)
	h = 1.0/(1 + n)^2
	f1 = (2x[1]-x[2]+0.5h^2*(x[1] + h + 1)^3)^2
	fn = (2x[n]-x[n-1]+0.5h^2*(x[n] + n*h + 1)^3)^2  
    return sum((2x[i]-x[i-1]-x[i+1]+0.5h^2*(x[i] + i*h + 1)^3)^2 for i in 2:n-1) + f1 + fn
end

g! = (g, x) -> begin
	r = similar(g)
	n = lastindex(x)
	h = 1.0/(1 + n)^2
	for i in 2:n-1
		r[i] = 2x[i]-x[i-1]-x[i+1]+0.5h^2*(x[i] + i*h + 1)^3
	end
	r[1] = 2x[1]-x[2]+0.5h^2*(x[1] + h + 1)^3
	r[n] = 2x[n]-x[n-1]+0.5h^2*(x[n] + n*h + 1)^3
	for i in 2:n-1
		g[i] = -r[i-1] -r[i+1] + r[i]*(2+3/2.0*h^2*(x[i] + i*h + 1)^2)
	end
	g[1] = -r[2] + r[1] *(2+3/2.0*h^2*(x[1] + h + 1)^2)
	g[n] = -r[n-1] + r[n] *(2+3/2.0*h^2*(x[n] + n*h + 1)^2)
    return g
end

fg! = (g, x) -> begin
	r = similar(g)
	n = lastindex(x)
	h = 1.0/(1 + n)^2
	fx = (2x[1]-x[2]+0.5h^2*(x[1] + h + 1)^3)^2 + (2x[n]-x[n-1]+0.5h^2*(x[n] + n*h + 1)^3)^2  
	for i in 2:n-1
		r[i] = 2x[i]-x[i-1]-x[i+1]+0.5h^2*(x[i] + i*h + 1)^3
		fx += r[i]^2
	end
	r[1] = 2x[1]-x[2]+0.5h^2*(x[1] + h + 1)^3
	r[n] = 2x[n]-x[n-1]+0.5h^2*(x[n] + n*h + 1)^3
	for i in 2:n-1
		g[i] = -r[i-1] -r[i+1] + r[i]*(2+3/2.0*h^2*(x[i] + i*h + 1)^2)
	end
	g[1] = -r[2] + r[1] *(2+3/2.0*h^2*(x[1] + h + 1)^2)
	g[n] = -r[n-1] + r[n] *(2+3/2.0*h^2*(x[n] + n*h + 1)^2)
	return fx, g
end

init = (n::Int=5000) -> begin
	n < 2 && @warn("MOREBV: number of variables must be ≥ 2")
	n = max(n, 2)
    return n, 0.5ones(n)
end
@warn "Debug MOREBV"
#TestSet["MOREBV"] = UncProgram("MOREBV",  f, g!, fg!, init)