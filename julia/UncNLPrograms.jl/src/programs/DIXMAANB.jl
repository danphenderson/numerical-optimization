#    Problem :
#    *********
#    The Dixon-Maany test problem (version B)
#
#    Origonal SIF Source:
#    L.C.W. Dixon and Z. Maany,
#    "A family of test problems with sparse Hessians for unconstrained
#    optimization",
#    TR 206, Numerical Optimization Centre, Hatfield Polytechnic, 1988.
#
# 	 See also Buckley#221 (p. 49)
#
#    DIXMAANB.SIF classification OUR2-AN-V-0
#
#    M is equal to the third of the number of variables (N is divisible by 3)
#
# Daniel Henderson, 08/2021

f = (x) -> begin
	n = lastindex(x)
	m = Int(n/3)
	fx = 1.0
	α = 1.0
	Β = 0.0625
	γ = 0.0625
	δ = 0.0625
	E = [0, 0, 0, 0]
	fx += sum(α*x[i]^2*(i/n)^E[1] for i in 1:n)							# S1
	fx += sum(Β*x[i]^2*(x[i+1]+x[i+1]^2)^2*(i/n)^E[2] for i in 1:n-1)   # S2
	fx += sum(γ*x[i]^2*x[i+m]^4*(i/n)^E[3] for i in 1:Int(2*m))   		# S3
	fx += sum(δ*x[i]*x[i+2m]*(i/n)^E[4] for i in 1:m) 					# S4
    return fx
end

g! = (g, x) -> begin
	n = lastindex(x)
	m = Int(n/3)
	α = 1.0
	Β = 0.0625
	γ = 0.0625
	δ = 0.0625
	E = [0, 0, 0, 0]

	# S1 := ∑ α*x[i]^2*(i/n)^E[1] for i in 1:n
	for i in 1:n
		g[i] += 2*α*x[i]*(i/n)^E[1]
	end

	# S2 := ∑ Β*x[i]^2*(x[i+1]+x[i+1]^2)^2*(i/n)^E[2] for i in 1:n-1
	for i in 1:n-1
		g[i]   += 2*Β*x[i]*(x[i+1]+x[i+1]^2)^2*(i/n)^E[2]
		g[i+1] += 2*Β*x[i]^2*(x[i+1]+x[i+1]^2)*(i/n)^E[2]*(1+2x[i+1])
	end
		
	# S3 := ∑ γ*x[i]^2*x[i+m]^4*(i/n)^E[3] for i in 1:Int(2*m)
	for i in 1:Int(2*m)
		g[i]   += 2*γ*x[i]*x[i+m]^4*(i/n)^E[3]
		g[i+m] += 4*γ*x[i]^2*x[i+m]^3*(i/n)^E[3]
	end

	#S4 := ∑ δ*x[i]*x[i+2m]*(i/n)^E[4] for i in 1:m
	for i in 1:m
		g[i]    += δ*x[i+2m]*(i/n)^E[4]
		g[i+2m] += δ*x[i]*(i/n)^E[4]
	end
	return g
end

fg! = (g, x) -> begin
	n = lastindex(x)
	m = Int(n/3)
	fx = 1.0
	α = 1.0
	Β = 0.0625
	γ = 0.0625
	δ = 0.0625
	E = [0, 0, 0, 0]

	# S1 := ∑ α*x[i]^2*(i/n)^E[1] for i in 1:n
	for i in 1:n
		g[i] += 2*α*x[i]*(i/n)^E[1]
		fx   += α*x[i]^2*(i/n)^E[1]
	end

	# S2 := ∑ Β*x[i]^2*(x[i+1]+x[i+1]^2)^2*(i/n)^E[2] for i in 1:n-1
	for i in 1:n-1
		g[i]   += 2*Β*x[i]*(x[i+1]+x[i+1]^2)^2*(i/n)^E[2]
		g[i+1] += 2*Β*x[i]^2*(x[i+1]+x[i+1]^2)*(i/n)^E[2]*(1+2x[i+1])
		fx     += Β*x[i]^2*(x[i+1]+x[i+1]^2)^2*(i/n)^E[2]
	end
		
	# S3 := ∑ γ*x[i]^2*x[i+m]^4*(i/n)^E[3] for i in 1:Int(2*m)
	for i in 1:Int(2*m)
		g[i]   += 2*γ*x[i]*x[i+m]^4*(i/n)^E[3]
		g[i+m] += 4*γ*x[i]^2*x[i+m]^3*(i/n)^E[3]
		fx     += γ*x[i]^2*x[i+m]^4*(i/n)^E[3]
	end

	#S4 := ∑ δ*x[i]*x[i+2m]*(i/n)^E[4] for i in 1:m
	for i in 1:m
		g[i]    += δ*x[i+2m]*(i/n)^E[4]
		g[i+2m] += δ*x[i]*(i/n)^E[4]
		fx      += δ*x[i]*x[i+2m]*(i/n)^E[4]
	end
	return fx, g
end

@warn "TODO: Clean up DIXMAANB (and others) implementation. See comment block"

init = (n::Int=3000) -> begin
	mod(n, 3) > 0 && @warn "DIXMAANB: number of variables must be divisible by 3" 
	q = max(1, div(n, 3))
	n = 3*q

	x0 = 2.0*ones(n)
    return n, x0
end

TestSet["DIXMAANB"] = UncProgram("DIXMAANB", f, g!, fg!, init)

# A cleaner implmentation might look like this: (issues rebasing S2, g[n] and g[m] need to be fixed)
# function g!(x, g)
# 	n = lastindex(x)
# 	m = Int(n/3)
# 	fx = 1.0
# 	α = 1.0
# 	Β = 0.0625
# 	γ = 0.0625
# 	δ = 0.0625
# 	E = [0, 0, 0, 0]
# 	for i in 1:m
# 		# S1 := ∑ α*x[i]^2*(i/n)^E[1] for i in 1:n
# 		g[i]    += 2*α*x[i]*(i/n)^E[1]
# 		g[i+m]  += 2*α*x[i]*(i/n)^E[1]
# 		g[i+2m] += 2*α*x[i]*(i/n)^E[1]

# 		# S2 := ∑ Β*x[i]^2*(x[i+1]+x[i+1]^2)^2*(i/n)^E[2] for i in 1:n-1
# 		g[i]    += 2*Β*x[i]*(x[i+1]+x[i+1]^2)^2*(i/n)^E[2]
# 		g[i+m]  += 2*Β*x[i]*(x[i+1]+x[i+1]^2)^2*(i/n)^E[2]
# 		g[i+2m] += 2*Β*x[i]*(x[i+1]+x[i+1]^2)^2*(i/n)^E[2]
# 		g[i]    += 2*Β*x[i-1]^2*(x[i]+x[i]^2)*((i-1)/n)^E[2]*(1+2x[i])
# 		g[i+m]  += 2*Β*x[i-1]^2*(x[i]+x[i]^2)*((i-1)/n)^E[2]*(1+2x[i])
# 		g[i+2m] += 2*Β*x[i-1]^2*(x[i]+x[i]^2)*((i-1)/n)^E[2]*(1+2x[i])
		
# 		# S3 := ∑ γ*x[i]^2*x[i+m]^4*(i/n)^E[3] for i in 1:Int(2*m)
# 		g[i]    += 2*γ*x[i]*x[i+m]^4*(i/n)^E[3]
# 		g[i+m]  += 2*γ*x[i]*x[i+m]^4*(i/n)^E[3]
# 		g[i+m]  += 4*γ*x[i]^2*x[i+m]^3*(i/n)^E[3]
# 		g[i+2m] += 4*γ*x[i]^2*x[i+m]^3*(i/n)^E[3]

# 		#S4 := ∑ δ*x[i]*x[i+2m]*(i/n)^E[4] for i in 1:m
# 		g[i]    += δ*x[i+2m]*(i/n)^E[4]
# 		g[i+2m] += δ*x[i]*(i/n)^E[4]
# 	end
# 	g[1] = 2*α*x[1]*n^(-E[1]) + 2*Β*x[1]*(x[2]+x[2]^2)^2*n^(-E[2]) + 2*γ*x[1]*x[1+m]^4*n^(-E[3]) + δ*x[2m+1]*n^(-E[4])
# 	g[n] -= 2*Β*x[n-1]^2*(x[n]+x[n]^2)*((n-1)/n)^E[2]*(1+2x[n]) 
# 	return g
# end