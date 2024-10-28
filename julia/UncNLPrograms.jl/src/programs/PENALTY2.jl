#    Problem : GROUP A
#    *********	 
#    The second penalty function
#
#    This is a nonlinear least-squares problem with M=2*N groups.
#    Group 1 is linear.
#    Groups 2 to N use 2 nonlinear elements.
#    Groups N+1 to M-1 use 1 nonlinear element.
#    Group M uses N nonlinear elements.
#    The Hessian matrix is dense.
#
#    Origonal SIF Source: Problem 24 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#
#    See also Buckley#112 (p. 80)
#
#    PENALTY2.SIF classification SUR2-AN-V-0
#
#    Number of variables is variable (n>3)
#
# Daniel Henderson, 08/2021   

f = (x) -> begin
	@warn "PENALTY2.jl doesn't agree with CUTEst or a JuMP (all three conflict)"
	N  = lastindex(x)
	fx = (x[1]-0.2)^2
	a  = 1e-5
	y  = ones(2*N)

  	for i in 1:2*N
    	y[i] = exp(i / 10.0) + exp((i - 1) / 10.0)
  	end

	α  = N*x[1]^2-1
	
	for i in 2:N
		fx += a*(exp(x[i]/10)+exp(x[i-1]/10) - y[i])^2
		fx += a*(exp(x[i]/10)-exp(-1/10))^2
		α += (N-i+1)*x[i]^2 - 1
	end
	
	fx += α^2
	return fx
end

# check chain rule use on product of sum term 
g! = (g, x) -> begin
    N  = lastindex(x)
	γ  = (1e-5)/5
	α  = N*x[1]^2-1

	g[1] = 2*(x[1]-0.2)*x[1]
	for i in 2:N
		g[i-1] += γ*exp(x[i-1]/10) * (exp(x[i]/10) + exp(x[i-1]/10) - exp((i+1)/10) - exp(i/10) )
		g[i]   -= γ*exp(x[i]/10) * (exp(1/10) + exp(i/10) + exp((i+1)/10) - exp(x[i-1]/10) - 2*exp(x[i]/10))
		α      += (N-i+1)*x[i]^2 - 1
	end
		α *= 2 # diff of outer function 
	for i in 2:N
		g[i] += α * 2*(N-i+1)*x[i] 
	end
	return g
end


fg! = (g, x) -> begin
    return error("PENALTY2 fg! not implmented")
end

init = (n::Int=1000) -> begin
	x0 = Vector(1.0:1.0:n)
    return n, x0
end

@warn "TODO: Debug PENALTY2"
#TestSet["PENALTY2"] = UncProgram("PENALTY2", f, g!, fg!, init)