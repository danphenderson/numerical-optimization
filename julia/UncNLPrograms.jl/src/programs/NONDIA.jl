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
#    NONDIA.jl.SIF classification OUR2-AN-V-0
#
#    Number of variables is variable
#
# Daniel Henderson, 08/2021  

f = (x) -> begin
	fx = (x[1]-1)^2
	for i in 2:lastindex(x)
		fx += 100*(x[1]-x[i-1]^2)^2
	end
    return fx
end

g! = (g, x) -> begin
	g[1] = 2*(x[1]-1)
	for i in 2:lastindex(x)
		g[1]   += 200(x[1]-x[i-1]^2)
		g[i-1] += -400x[i-1]*(x[1] - x[i-1]^2)
	end
    return g
end

fg! = (g, x) -> begin
	fx = (x[1]-1)^2
	g[1] = 2*(x[1]-1)
	for i in 2:lastindex(x)
		fx += 100*(x[1]-x[i-1]^2)^2
		g[1]   += 200(x[1]-x[i-1]^2)
		g[i-1] += -400x[i-1]*(x[1] - x[i-1]^2)
	end
    return g
end

init = (n::Int=5000) -> begin
	n < 2 && @warn("ARWHEAD: number of variables must be ≥ 2")
	n = max(n, 2)
    
    return n, -1.0ones(n)
end

TestSet["NONDIA"] = UncProgram("NONDIA", f, g!, fg!, init)