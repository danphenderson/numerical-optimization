#    Problem : GROUP A
#    *********
#	 A nondiagonal quartic test problem.
#	 
#	 This problem has an arrow-head type Hessian with a tridiagonal
#    central part and a border of width 1.
#    The Hessian is singular at the solution.
#
#    Origonal SIF Source:
#    A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
#    "Performance of a multi-frontal scheme for partially separable
#    optimization"
#    Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
#
#    NONDQUAR.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable
#
# Daniel Henderson, 08/2021  

f = (x) -> begin
	n  = lastindex(x)
	fx =  (x[1]-x[2])^2 + (x[n-1] - x[n])^2
	for i in 1:n-2
      	fx += (x[i]+x[i+1]+x[n])^4
	end
    return fx
end

g! = (g, x) -> begin
	n    = lastindex(x)
	g[1] = 2(x[1] - x[2]) + 4(x[1] + x[2] + x[n])^3
	g[2] = -2(x[1] - x[2]) + 4(x[1] + x[2] + x[n])^3 + 4(x[2] + x[3] + x[n])^3
	gn   = 4(x[1] + x[2] + x[n])^3 +  4(x[2] + x[3] + x[n])^3
	for i in 3:n-2
		g[i] += 4(x[i] + x[i+1] + x[n])^3 + 4(x[i+1] + x[i+2] + x[n])^3
		gn   += 4(x[i] + x[i+1] + x[n])^3
	end
	g[n-1] = 2(x[n-1] - x[n]) + 4(x[n-2] + x[n-1] + x[n])^3
	g[n]   = gn - 2(x[n-1] - x[n])
	return g
end

fg! = (g, x) -> begin
	n    = lastindex(x)
	fx   =  (x[1]-x[2])^2 + (x[n-1] - x[n])^2 + (x[1]+x[2]+x[n])^4 + (x[2]+x[3]+x[n])^4
	g[1] = 2(x[1] - x[2]) + 4(x[1] + x[2] + x[n])^3 
	g[2] = -2(x[1] - x[2]) + 4(x[1] + x[2] + x[n])^3 + 4(x[2] + x[3] + x[n])^3
	gn   = 4(x[1] + x[2] + x[n])^3 +  4(x[2] + x[3] + x[n])^3
	for i in 3:n-2
		g[i] += 4(x[i] + x[i+1] + x[n])^3 + 4(x[i+1] + x[i+2] + x[n])^3
		gn   += 4(x[i] + x[i+1] + x[n])^3
		fx += (x[i]+x[i+1]+x[n])^4
	end
	g[n-1] = 2(x[n-1] - x[n]) + 4(x[n-2] + x[n-1] + x[n])^3
	g[n]   = gn - 2(x[n-1] - x[n])
	return fx, g
end

init = (n::Int=5000) -> begin
	mod(n, 2) > 0 && @warn "NONDQUAR: number of variables must be even" 
	q = max(1, div(n, 2))
	n = 2q
	x0 = [j % 2 == 0 ? -1.0 : 1.0 for j in 1:n]
    return n, x0
end

#TestSet["NONDQUAR"] = UncProgram("NONDQUAR", f, g!, fg!, init)
