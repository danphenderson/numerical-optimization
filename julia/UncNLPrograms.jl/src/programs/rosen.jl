f = x -> begin
	N = lastindex(x)
    return 100sum((x[i+1] - x[i]^2)^2 for i = 1:N-1) + sum((x[i] - 1.0)^2 for i = 1:N-1)
end

g! = (g, x) -> begin
	N = lastindex(x)
	g[1] = -2*(1 - x[1]) - 400*x[1]*(-x[1]^2 + x[2])
	for i in 2:N-1
		g[i] = -2*(1 - x[i]) + 200*(-x[i-1]^2 + x[i]) - 400*x[i]*(-x[i]^2 + x[1 + i])
	end
	g[N] = 200 * (x[N] - x[N-1]^2)
	return g
end

fg! = (g, x) -> begin
	@warn "rosen fg! not implmented"
end

init = (n::Int=500) -> begin
	x0 = [j/(n+1) for j in 1:n]
	return n, x0
end

TestSet["rosen"] = UncProgram("rosen",  f, g!, fg!, init)
