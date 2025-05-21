using LinearAlgebra
"""
Implementing Algorithm 16.5 From Section 16.7
in Nocedal and Wright's text on Numerical Optimization
"""

mutable struct Constraints
    active_set::Set
    # TODO
end

function _is_feasible(q, xk, constraints)
    # TODO
end

function _is_kkt_satisfied(q::Function, xk::Vector, constraints::Constraints)
    # TODO
end

function _find_cauchy_point(q::Function, xk::Vector, constraints::Constraints)
    # 
end

function _solve_subproblem(q::Function, xk::Vector, constraints::Constraints)
    # TODO here we confirm in leads to decrease in `q` and that `x^+` is feasible
end

function projected_gradients(q::Function, x0::Vector, constraints::Constraints)
    xk = x0
    for i ∈ range(1, 1000)
        if _is_kkt_satisfied(q, x0) # terminal condition
            return x0
        end # otherwise solve subproblem
        xc = _find_cauchy_point(q, xk, constraints)
        xk = _solve_subproblem(q, xk, constraints)
    end
end