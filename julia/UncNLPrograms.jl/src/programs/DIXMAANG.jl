#    Problem :
#    *********
#    The Dixon-Maany test problem (version G)
#
#    Origonal SIF Source:
#    L.C.W. Dixon and Z. Maany,
#    "A family of test problems with sparse Hessians for unconstrained
#    optimization",
#    TR 206, Numerical Optimization Centre, Hatfield Polytechnic, 1988.
#
# 	 See also Buckley#221 (p. 49)
#
#    DIXMAANG.SIF classification OUR2-AN-V-0
#
#    m is equal to the third of the number of variables (n is divisible by 3)
#
# Daniel Henderson, 08/2021

f = (x) -> begin
	n = lastindex(x)
	m = Int(n/3)
	fx = 1.0

	α = 1.0
	Β = 0.125
	γ = 0.125
	δ = 0.125
	E = [1, 0, 0, 1]
	
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
	Β = 0.125
	γ = 0.125
	δ = 0.125
	E = [1, 0, 0, 1]

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
	Β = 0.125
	γ = 0.125
	δ = 0.125
	E = [1, 0, 0, 1]

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

init = (n::Int=3000) -> begin
	mod(n, 3) > 0 && @warn "DIXMAANG: number of variables must be divisible by 3" 
	q = max(1, div(n, 3))
	n = 3*q

	x0 = 2.0*ones(n)
    return n, x0
end

TestSet["DIXMAANG"] = UncProgram("DIXMAANG", f, g!, fg!, init)