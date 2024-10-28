#    Problem : GROUP A
#    *********
#	 The extended Woods problem.
#	 This problem is a sum of n/4 sets of 6 terms, each of which is
#    assigned its own group.  For a given set i, the groups are
#    A(i), B(i), C(i), D(i), E(i) and F(i). Groups A(i) and C(i) contain 1
#    nonlinear element each, denoted Y(i) and Z(i).
#
#    The problem dimension is defined from the number of these sets.
#    The number of problem variables is then 4 times larger.
#
#	 This version uses a slightly unorthodox expression of Woods
#    function as a sum of squares (see Buckley)
#
#    Origonal SIF Source: problem 14 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#
#    See also Toint#27, Buckley#17 (p. 101), Conn, Gould, Toint#7
#
#    WOODS.SIF classification SUR2-AN-V-0
#
#	 This problem is decomposed in n linear groups, the last n-1 of which
#    are 2 x 2 and singular.
#
#    NS is the number of sets (= n/4)
#
# Daniel Henderson, 08/2021  

f = (x) -> begin
	fx = 0.0
	for i in 1:Int(lastindex(x)/4)
		fx += (100*(x[4*i-2]-x[4*i-3]^2)^2 + (1-x[4*i-3])^2 + 90*(x[4*i]-x[4*i-1]^2)^2 + (1-x[4*i-1])^2 + 10*(x[4*i-2]+x[4*i]-2)^2 + 0.1*(x[4*i-2]-x[4*i])^2)
	end
    return fx
end

g! = (g, x) -> begin
	for i in 1:Int(lastindex(x)/4)
		g[4i-3] = -2(1-x[4i-3]) - 400*x[4i-3]*(x[4i-2]-x[4i-3]^2)
		g[4i-2] = 200(x[4i-2]-x[4i-3]^2) + 0.2(x[4i-2]-x[4i]) + 20(x[4i-2] + x[4i]-2)
		g[4i-1] = -2(1-x[4i-1]) - 360*x[4i-1]*(x[4i]-x[4i-1]^2)
		g[4i] = -0.2(x[4i-2]-x[4i]) + 20(x[4i-2] + x[4i]-2) + 180(x[4i]-x[4i-1]^2)
	end
    return g
end

fg! = (g, x) -> begin
	fx = 0.0
	for i in 1:Int(lastindex(x)/4)
		fx += (100*(x[4*i-2]-x[4*i-3]^2)^2 + (1-x[4*i-3])^2 + 90*(x[4*i]-x[4*i-1]^2)^2 + (1-x[4*i-1])^2 + 10*(x[4*i-2]+x[4*i]-2)^2 + 0.1*(x[4*i-2]-x[4*i])^2)
		g[4i-3] = -2(1-x[4i-3]) - 400*x[4i-3]*(x[4i-2]-x[4i-3]^2)
		g[4i-2] = 200(x[4i-2]-x[4i-3]^2) + 0.2(x[4i-2]-x[4i]) + 20(x[4i-2] + x[4i]-2)
		g[4i-1] = -2(1-x[4i-1]) - 360*x[4i-1]*(x[4i]-x[4i-1]^2)
		g[4i] = -0.2(x[4i-2]-x[4i]) + 20(x[4i-2] + x[4i]-2) + 180(x[4i]-x[4i-1]^2)
	end
    return fx, g
end

init = (n::Int=4000) -> begin
	mod(n, 4) > 0 && @warn "WOODS: number of variables must be divisible by 4" 
	q = max(1, div(n, 4))
	n = 4*q

	x0 = [j % 2 == 1 ? -3.0 : -1.0 for j in 1:n]
    return n, x0
end

TestSet["WOODS"] = UncProgram("WOODS", f, g!, fg!, init)