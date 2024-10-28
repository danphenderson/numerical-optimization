#    Problem :
#    *********	 
#
#    The Boundary Value problem.
#    This is the nonlinear least-squares version without fixed variables.
#  
#
#    Source:  problem 28 in
#    J.J. More', B.S. Garbow and K.E. Hillstrom,
#    "Testing Unconstrained Optimization Software",
#    ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
#	
#    Implementation translated from:
#    See also Buckley#17 (p. 94):
#	
#	
#    SIF input: Ph. Toint, Dec 1989 and Nick Gould, Oct 1992.
#
#    classification SUR2-MN-V-0
#    The number of variables is N. 
#
# Daniel Henderson, 08/2021   

# f = x -> begin

# end

# g! = (g, x) -> begin
# 	n = lastindex(x)
# end

# fg! = (g, x) -> begin

# end

# init = (n::Int=5000) -> begin

# end

# @warn "FRUEROTH needs to be implemented"

# TestSet["FREUROTH"] = UncProgram("FREUROTH",  f, g!, fg!, init)