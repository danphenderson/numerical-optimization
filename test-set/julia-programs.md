```julia
f = x -> begin
    N = lastindex(x)
    return sum((x[i] - 1.0)^2 for i = 1:N) + (sum(x[j]^2 for j = 1:N) - 0.25)^2
end
init = (n::Int=500) -> begin
    x0 = [j for j in 1:n]
    return n, x0, "Extended Penalty"
end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(i * x[i]^2 for i = 1:N) + (1/100) * (sum(x[i] for i = 1:N))^2
 end
 init = (n::Int=500) -> begin
     x0 = [0.5 for _ in 1:n]
     return n, x0, "Perturbed Quadratic"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(i/10 * (exp(x[i]) - x[i]) for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Raydan 1"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(exp(x[i]) - x[i] for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Raydan 2"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(exp(x[i]) - i * sin(x[i]) for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Diagonal 3"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(exp(x[i]) - i * x[i] for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Hager"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(x[i]^2 + x[i] * x[i+1] - 4 * x[i+1]^2 + 3 * x[i+1] for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = [2.0 for _ in 1:n]
     return n, x0, "Generalized Tridiagonal 1"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(5*(x[i] - x[i+1])^2 + 3*(x[i] - x[i+1])^2 + 11*(x[i] - x[i+1]) for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Generalized Tridiagonal 2"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(log(exp(x[i]) + exp(-x[i])) for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = [1.1 for _ in 1:n]
     return n, x0, "Diagonal 5"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return 100sum((x[i+1] - x[i]^3)^2 for i = 1:N-1) + sum((1.0 - x[i]^2)^2 for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = [(-1.2 + 2*(mod(j, 2))) for j in 1:n]
     return n, x0, "Generalized White & Holst"
 end
```
```julia
 f = x -> begin
     N = lastindex(x) - 1
     return sum((x[i]^2 + x[i+1]^2 + x[i]*x[i+1])^2 + sin(x[i])^2 + cos(x[i])^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = [ifelse(mod(j, 2) == 1, 3.0, 0.1) for j in 1:n]
     return n, x0, "Generalized PSC1"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (x[1] - 3.0)^2 + sum((x[1] - 3.0 - 2 * sum(x[j] for j = 1:i)^2)^2 for i = 2:N)
 end
 init = (n::Int=500) -> begin
     x0 = [0.01 for _ in 1:n]
     return n, x0, "Full Hessian FH1"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (x[1] - 5)^2 + sum((sum(x[1:i]) - 1)^2 for i = 2:N)
 end
 init = (n::Int=500) -> begin
     x0 = [0.01 for _ in 1:n]
     return n, x0, "Full Hessian FH2"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (sum(x[i] for i = 1:N))^2 + sum((i/100.0) * x[i]^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = [0.5 for _ in 1:n]
     return n, x0, "Perturbed Quadratic Diagonal"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return 0.5 * sum(i * x[i]^2 for i = 1:N) - x[N]
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Quadratic QF1"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum((x[i]^2 - 2.0)^2 for i = 1:N-1) + sum((x[i]^2 - 0.5)^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Extended quadratic penalty QP1"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (sum(x[i]^2 for i = 1:N) - 100)^2 + sum((x[i]^2 - sin(x[i]))^2 for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Extended quadratic penalty QP2"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return 0.5 * sum(i * (x[i]^2 - 1)^2 for i = 1:N) - x[N]
 end
 init = (n::Int=500) -> begin
     x0 = [0.5 for _ in 1:n]
     return n, x0, "Quadratic QF2"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum((x[i]*x[i+1] - 1.0)^2 + 0.1*(x[i] + 1.0)*(x[i+1] + 1.0) for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Extended Tridiagonal 2"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     p = 1 / 10^8
     h = 1 / (N + 1)
     c = 1
     term1 = 0.5 * p * (x[1]^2 + x[N]^2)
     term2 = sum((x[i] - x[i+1])^2 for i = 1:N-1)
     term3 = sum((p * (h^2 + 2) / h^2) * x[i] + (c * p / h^2) * cos(x[i]) for i = 1:N)
     return term1 + term2 - term3
 end
 init = (n::Int=500) -> begin
     h = 1 / (n + 1)
     x0 = [j * h for j in 1:n]
     return n, x0, "FLETCBV3 (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return 100 * sum((x[i+1] - x[i] + 1 - x[i]^2)^2 for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = zeros(n)
     return n, x0, "FLETCHCR (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     term1 = sum((-4 * x[i] + 3)^2 for i = 1:N-4)
     term2 = (x[1]^2 + 2 * x[2]^2 + 3 * x[3]^2 + 4 * x[4]^2 + 5 * x[N]^2)^2
     return term1 + term2
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "BDQRTIC (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     alpha, beta, gamma, delta = 2, 1, 1, 1
     return gamma * (delta * x[1] - 1)^2 + sum(i * (alpha * x[i] - beta * x[i-1])^2 for i = 2:N)
 end
 init = (n::Int=500) -> begin
     x0 = [1.0 for _ in 1:n]
     return n, x0, "TRIDIA (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(i * j * x[j] - 1 for i = 1:N, j = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "ARGLINB (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(-4*x[i] + 3 for i = 1:N-1) + sum((x[i]^2 + x[N]^2)^2 for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "ARWHEAD (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (x[1] - 1.0)^2 + sum(100 * (x[1] - x[i-1]^2)^2 for i = 2:N)
 end
 init = (n::Int=500) -> begin
     x0 = [-1.0 for _ in 1:n]
     return n, x0, "NONDIA (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (x[1] - x[2])^2 + sum((x[i] + x[i+1] + x[N])^4 for i = 1:N-2) + (x[N-1] + x[N])^2
 end
 init = (n::Int=500) -> begin
     x0 = [(-1)^j for j in 1:n]
     return n, x0, "NONDQUAR (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     c, d = 100, 100
     return sum(x[i]^2 + c * x[i+1] + d * x[i+2]^2 for i = 1:N-2)
 end
 init = (n::Int=500) -> begin
     x0 = [3.0 for _ in 1:n]
     return n, x0, "DQDRTIC (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(sin(x[1] + x[i]^2 - 1) for i = 1:N-1) + 0.5 * sin(x[N]^2)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "EG2 function (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     k = 20
     return sum((q_i(x, i, N, k)^4 - 20q_i(x, i, N, k) - 0.1q_i(x, i, N, k)) for i = 1:N)
 end
 q_i = (x, i, N, k) -> begin
     if i <= N - k
         return sum(x[j] for j = i:i+k)
     else
         return sum(x[j] for j = i:N)
     end
 end
 init = (n::Int=500) -> begin
     x0 = [0.001 / (n+1) for _ in 1:n]
     return n, x0, "CURLY20"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return x[1]^2 + sum(i * x[i]^2 + (1/100) * sum(x[j] for j = 1:i)^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = [0.5 for _ in 1:n]
     return n, x0, "Partially Perturbed Quadratic"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     term1 = (3x[1] - 2x[1]^2)^2
     termN = (3x[N] - 2x[N]^2 - x[N-1] + 1)^2
     return term1 + sum((3x[i] - 2x[i]^2 - x[i-1] - 2x[i+1] + 1)^2 for i = 2:N-1) + termN
 end
 init = (n::Int=500) -> begin
     x0 = [-1.0 for _ in 1:n]
     return n, x0, "Broyden Tridiagonal"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(i * x[i]^2 for i = 1:N) + (1/100) * (x[1] + x[N])^2
 end
 init = (n::Int=500) -> begin
     x0 = fill(0.5, n)
     return n, x0, "Almost Perturbed Quadratic"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return x[1]^2 + sum(i * x[i]^2 - (x[i-1] + x[i] + x[i+1])^2 for i = 2:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = [0.5 for _ in 1:n]
     return n, x0, "Perturbed Tridiagonal Quadratic"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum((sum(x[j] for j = 1:i))^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Staircase 1"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum((sum(x[j] for j = 1:i) - i)^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = zeros(n)
     return n, x0, "Staircase 2"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return 4sum((-x[1] + x[i]^2)^2 for i = 1:N) + sum((x[i] - 1.0)^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = [4.0 for _ in 1:n]
     return n, x0, "LIARWHD (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum((i * x[i])^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "POWER (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum((x[i]^2 + x[i+1]^2)^2 for i = 1:N-1) + sum(-4 * x[i] + 3 for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = [2.0 for _ in 1:n]
     return n, x0, "ENGVAL1"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return 16 + sum((x[i] - 2)^4 + (x[i] * x[i+1] - 2 * x[i+1])^2 + (x[i+1] + 1)^2 for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = zeros(n)
     return n, x0, "EDENSCH (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(cos(2*x[i] - x[N] - x[1]) for i = 2:N-1) + sum(x[i] for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = [j/(n+1) for j in 1:n]
     return n, x0, "INDEF (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (x[1] - 1)^2 + sum(100 * (x[i] - x[i-1]^3)^2 for i = 2:N)
 end
 init = (n::Int=500) -> begin
     x0 = [(-1.2)^mod(j, 2) * 1.0 for j in 1:n]  # Alternates between -1.2 and 1
     return n, x0, "CUBE (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(exp(0.1 * x[i] * x[i+1]) for i = 1:N-1) - 101 * sum(i * x[i] for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = fill(0.0, n)
     return n, x0, "EXPLIN1 (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum((x[i] + x[i+1]) * exp(-x[i+2] * (x[i] + x[i+1])) for i = 1:N-2)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "BDEXP (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     term1 = (sum(x[i] for i = 1:N))^2
     term2 = sum(x[i] + 0.5 * x[i]^2 for i = 1:N)
     term3 = 2 * sum((sum(x[j] for j = i:N))^2 for i = 2:N)
     return term1 - term2 + term3
 end
 init = (n::Int=500) -> begin
     x0 = [j for j in 1:n]
     return n, x0, "HARKERP2"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(sin(2x[i])^2 * sin(2x[i+1])^2 + 0.05 * (x[i]^2 + x[i+1]^2)^2 for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = [(-1)^j * 506.2 for j in 1:n]
     return n, x0, "GENHUMPS (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(-1.5 * x[i] + 2.5 * x[i+1]^2 + 1 + (x[i] - x[i+1])^2 + sin(x[i] + x[i+1]) for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = [1.0 for j in 1:n]
     return n, x0, "MCCORMCK (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (x[1] - 1.0)^2 + sum(4 * (x[i] - x[i-1]^2)^2 for i = 2:N)
 end
 init = (n::Int=500) -> begin
     x0 = [3.0 for _ in 1:n]
     return n, x0, "NONSCOMP (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum((x[i] - 1.0)^2 * (i / N) for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = [j for j in 1:n]
     return n, x0, "VARDIM (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum((x[i] - 1.0)^4 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = [2.0 for _ in 1:n]
     return n, x0, "QUARTC"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(exp(-x[i]) - (1 - x[i]) for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Diagonal 6"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     term1 = (x[1] - 1)^4
     term2 = sum((sin(x[i] - x[N]) - x[1]^2 + x[i]^2)^2 for i = 2:N-1)
     term3 = (x[N]^2 - x[1]^2)^2
     return term1 + term2 + term3
 end
 init = (n::Int=500) -> begin
     x0 = [0.1 for _ in 1:n]
     return n, x0, "SINQUAD (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return 4sum((x[i]^2 - x[1])^2 for i = 1:N) + sum((x[i] - 1.0)^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = [4.0 for _ in 1:n]
     return n, x0, "LIARWHD (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (x[1] - 1.0)^2 + sum((x[i] - x[i+1])^2 for i = 1:N-1) + (x[N] - 1.0)^2
 end
 init = (n::Int=500) -> begin
     x0 = [-1.0 for _ in 1:n]
     return n, x0, "DIXON3DQ"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(cos(-0.5 * x[i+1] + x[i]^2) for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "COSINE (CUTE)"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(sin(-0.5 * x[i+1] + x[i]^2) for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "SINE"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (x[1] - 1.0)^2 + sum((x[i+1] - x[i])^2 for i = 1:N-1) + (1.0 - x[N])^2
 end
 init = (n::Int=500) -> begin
     x0 = zeros(n)
     return n, x0, "BIGGSB1"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(x[i]^2 + (x[i+1] + x[i]^2)^2 for i = 1:N-1)
 end
 init = (n::Int=500) -> begin
     x0 = [1.0 for _ in 1:n]
     return n, x0, "Generalized Quartic"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(exp(x[i]) - 2x[i] - x[i]^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Diagonal 7"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return sum(x[i] * exp(x[i]) - 2 * x[i] - x[i]^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Diagonal 8"
 end
```
```julia
 f = x -> begin
     N = lastindex(x)
     return (sum(x[i] for i = 1:N))^2 + sum(x[i]*exp(x[i]) - 2*x[i] - x[i]^2 for i = 1:N)
 end
 init = (n::Int=500) -> begin
     x0 = ones(n)
     return n, x0, "Full Hessian FH3"
 end
```
