#    Problem : GROUP A
#    *********
#    Source: problem 10 in
#    Ph.L. Toint,
#    "Test problems for partially separable optimization and results
#    for the routine PSPMIN",
#    Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
#
#    See also Buckley#116.  Note that MGH#21 is the separable version.
#    SIF input: Ph. Toint, Dec 1989.
#
#    problem 29 in
#    L. Luksan, C. Matonoha and J. Vlcek
#    Modified CUTE problems for sparse unconstrained optimization,
#    Technical Report 1081,
#    Institute of Computer Science,
#    Academy of Science of the Czech Republic
#
#    http://www.cs.cas.cz/matonoha/download/V1081.pdf
#
#    classification SUR2-AN-V-0
#
# # Daniel Henderson, 09/2021 


f = (x) -> begin # this is correct
    return 100.0*sum((x[i] - x[i - 1]^2)^2 for i = 1:2:(lastindex(x)-1)) + (1.0 - x[1])^2
end


g! = (g, x) -> begin
    @warn "EXTROSEN g! not implemented"
end


fg! = (g, x) -> begin
    @warn "EXTROSEN g! not implemented"
end

init = (n::Int=1000) -> begin
	n < 2 && @warn("EXTROSNB: number of variables must be ≥ 2")
	n = max(n, 2)
    
    return n, -1.0*ones(n)
end

#TestSet["EXTROSNB"] = UncProgram("EXTROSNB",  f, g!, fg!, init)
