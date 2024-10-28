#   Problem :
#   *********
#
#   Source: Implementation found here
#   https://www.sfu.ca/~ssurjano/powell.html
#
#   N is the number of free variables, and should be a multiple of 4
#
# Daniel Henderson, 08/2021

f = (x) -> begin  
    return error("LANGER.jl not implemented")
end

g! = (g, x) -> begin  
    return error("LANGER.jl not implemented")
end

fg! = (g, x) -> begin  
    return error("LANGER.jl not implemented")
end

@warn "TODO: Figure out why LANGER doesn't have an SIF but in Princeton and Georgian repo"

init = (n::Int=3000) -> begin 
    x0 = zeros(n)
    return n, x0
end

# TestSet["LANGER"] = UncProgram("LANGER", f, g!, fg!, 25, zeros(25))