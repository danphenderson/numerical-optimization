#    Problem : 
#    *********	 
#   The generalized Rosenbrock function.
#
#   Source: problem 5 in
#   S. Nash,
#   "Newton-type minimization via the Lanczos process",
#
#   SIF input: Nick Gould, Oct 1992.
#              minor correction by Ph. Shott, Jan 1995.
#
#   classification SUR2-AN-V-0
#
#   Number of variables
#
#
# Daniel Henderson, 08/2021   



f = x -> begin
	N = lastindex(x)
    return 1.0 + 100sum((x[i] - x[i-1]^2)^2 for i = 2:N) + sum((x[i] - 1.0)^2 for i = 2:N)
end

g! = (g, x) -> begin
	N = lastindex(x)
	g[1] = -400 * x[1] * (x[2] - x[1]^2)
	α    =  200 * (x[2] - x[1]^2)
	for i in 2:N-1
		g[i] =  -400 * x[i] * (x[i+1] - x[i]^2) + 2 * (x[i] - 1.0) + α
		α    =   200 * (x[i+1] - x[i]^2)
	end
	g[N] = 200 * (x[N] - x[N-1]^2) + 2*(x[N] - 1.0)
	return g
end

fg! = (g, x) -> begin
	@warn "GENROSE fg! not implmented"
	# N = lastindex(x)
	# fx = 0.0
	# g[1] = -400 * x[1] * (x[2] - x[1]^2)
	# α    =  200 * (x[2] - x[1]^2)
	# for i in 2:N-1
	# 	g[i] =  -400 * x[i] * (x[i+1] - x[i]^2) + 2 * (x[i] - 1.0) + α
	# 	α    =   200 * (x[i+1] - x[i]^2)
	# 	fx  += (x[i] - x[i-1]^2)^2 + (x[i] - 1.0)^2
	# end
	# fx  += (x[N] - x[N-1]^2)^2 + (x[N] - 1.0)^2
	# g[N] = 200 * (x[N] - x[N-1]^2) + 2*(x[N] - 1.0)
	# return fx,  g	
end

init = (n::Int=500) -> begin
	x0 = [j/(n+1) for j in 1:n]
	return n, x0
end

TestSet["GENROSE"] = UncProgram("GENROSE",  f, g!, fg!, init)