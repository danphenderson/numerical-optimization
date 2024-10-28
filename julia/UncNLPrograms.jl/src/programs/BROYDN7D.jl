#    Problem :
#    *********
#    The BROYDN7D problem by Oren.
#
#    Source:
#    Ph.L. Toint,
#    "Some numerical results using a sparse matrix updating formula in
#    unconstrained optimization",
#    Mathematics of Computation, vol. 32(114), pp. 839-852, 1978.
#
#    See also Buckley#84
#
#    BROYDN7D.SIF classification OUR2-AN-V-0
#
#    Number of variables must be greater than 4
#
#    This problem has a reputation of issues.
#
# Daniel Henderson, 08/2021  



f = (x) -> begin
	n = length(x)
	fx = abs(1-2x[2]+(3-2x[1])x[1])^(7/3) + abs(1-x[n-1]+(3-2x[n])x[n])^(7/3) 
	for i in 1:Int(n/2)
		fx += abs(x[i]+x[i+Int(n/2)])^(7/3)
	end
	for i in 2:n-1
		fx += abs(1-x[i-1]-2*x[i+1]+(3-2*x[i])*x[i])^(7/3)
	end
	return fx
end

g! = (g, x) -> begin
	n      = length(x)
	first  = -2.0*x[1] + 1 + (3.0-2.0*x[1])x[1]
	last   = -x[n-1] + 1 + (3.0-2.0*x[n])x[n]
	g[1]   = 7/3*abs(first)^(4/3) * ( first > 0 ? 1 : -1) * (3.0-4x[1])
	g[2]   = -14/3*abs(first)^(4/3) * ( first > 0 ? 1 : -1)
	g[n-1] = -7/3*abs(last)^(4/3) * ( last > 0 ? 1 : -1)
	g[n]   = 7/3*abs(last)^(4/3) * ( last > 0 ? 1 : -1) * (3.0-4.0x[n])
	last   = x[1] + x[Int(n/2)+1]
	g[1]  += 7/3*abs(last)^(4/3) * ( last > 0 ? 1 : -1) 
	g[Int(n/2)+1]  += 7/3*abs(last)^(4/3) * ( last > 0 ? 1 : -1)

	for i in 2:Int(n/2)
		first   = 1 - x[i-1] - 2x[i+1] + (3-2x[i])x[i]
		last    = x[i] + x[i + Int(n/2)]
		g[i-1] += -7.0/3*abs(first)^(4/3) * (first > 0 ? 1 : -1)
		g[i]   += 7.0/3*abs(first)^(4/3) * (first > 0 ? 1 : -1)*(3.0-4.0x[i]) + 7.0/3*abs(last)^(4/3) * (last > 0 ? 1 : -1)
		g[i+1] += -14/3*abs(first)^(4/3) * (first > 0 ? 1 : -1 )
		g[i+Int(n/2)] += 7.0/3*abs(last)^(4/3) * (last > 0 ? 1 : -1)
	end

	for i in Int(n/2)+1:n-1
		first   = 1 - x[i-1] - 2x[i+1] + (3-2x[i])x[i]
		g[i-1] += -7.0/3*abs(first)^(4/3) * (first > 0 ? 1 : -1)
		g[i]   += 7.0/3*abs(first)^(4/3) * (first > 0 ? 1 : -1)*(3.0-4.0x[i])
		g[i+1] += -14/3*abs(first)^(4/3) * (first > 0 ? 1 : -1 )
	end

	return g
end

fg! = (g, x) -> begin
	n      = length(x)
	first  = -2.0*x[1] + 1 + (3.0-2.0*x[1])x[1]
	last   = -x[n-1] + 1 + (3.0-2.0*x[n])x[n]
	fx     = abs(first)^(7/3) + abs(last)^(7/3)
	g[1]   = 7/3*abs(first)^(4/3) * ( first > 0 ? 1 : -1) * (3.0-4x[1])
	g[2]   = -14/3*abs(first)^(4/3) * ( first > 0 ? 1 : -1)
	g[n-1] = -7/3*abs(last)^(4/3) * ( last > 0 ? 1 : -1)
	g[n]   = 7/3*abs(last)^(4/3) * ( last > 0 ? 1 : -1) * (3.0-4.0x[n])
	last   = x[1] + x[Int(n/2)+1]
	fx    += abs(last)^(7/3)
	g[1]  += 7/3*abs(last)^(4/3) * ( last > 0 ? 1 : -1) 
	g[Int(n/2)+1]  += 7/3*abs(last)^(4/3) * ( last > 0 ? 1 : -1)

	for i in 2:Int(n/2)
		first  = 1 - x[i-1] - 2x[i+1] + (3-2x[i])x[i]
		last   = x[i] + x[i + Int(n/2)]
		fx 	   += abs(first)^(7/3) + abs(last)^(7/3)
		g[i-1] += -7.0/3*abs(first)^(4/3) * (first > 0 ? 1 : -1)
		g[i]   += 7.0/3*abs(first)^(4/3) * (first > 0 ? 1 : -1)*(3.0-4.0x[i]) + 7.0/3*abs(last)^(4/3) * (last > 0 ? 1 : -1)
		g[i+1] += -14/3*abs(first)^(4/3) * (first > 0 ? 1 : -1 )
		g[i+Int(n/2)] += 7.0/3*abs(last)^(4/3) * (last > 0 ? 1 : -1)
	end

	for i in Int(n/2)+1:n-1
		first   = 1 - x[i-1] - 2x[i+1] + (3-2x[i])x[i]
		fx 	   += abs(first)^(7/3)
		g[i-1] += -7.0/3*abs(first)^(4/3) * (first > 0 ? 1 : -1)
		g[i]   += 7.0/3*abs(first)^(4/3) * (first > 0 ? 1 : -1)*(3.0-4.0x[i])
		g[i+1] += -14/3*abs(first)^(4/3) * (first > 0 ? 1 : -1 )
	end

	return fx, g
end

init = (n::Int=5000) -> begin
	mod(n, 2) > 0 && @warn "BROYDN7D: number of variables must be even" 
	q = max(1, div(n, 2))
	n = 2q

    return n, ones(n)
end

@warn "TODO: Issue in CUTEst? See comment https://github.com/JuliaSmoothOptimizers/OptimizationProblems.jl/blob/main/src/broydn7d.jl#L50 "
# TestSet["BROYDN7D"] = UncProgram("BDQRTIC", f, g!, fg!, init)