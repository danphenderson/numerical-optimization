using CSV, DataFrames, OptimizationProblems, ADNLPModels, NLPModels
using LinearAlgebra, Statistics, Distributions
using OptimizationProblems
using JSON
using DataFrames
using CSV
using Random


Random.seed!(1234)

get_test_set_names() = 	[
    "genrose",
    "arglina",
    "freuroth",
    "eg2",
    "cosine",
    "arglinb",
    "arglinc",
    "argtrig",
    "arwhead",
    "bdqrtic",
    "brownal",
    "broyden3d",
    "chnrosnb_mod",
    "cragglvy",
    "cragglvy2",
    "curly10",
    "curly10",
    "curly20",
    "curly30",
    "dixon3dq",
    "dqdrtic",
    "dqrtic",
    "edensch",
    "engval1",
    "errinros_mod",
    "extrosnb",
    "fletcbv2",
    "fletcbv3_mod",
    "fletchcr",
    "genhumps",
    "genrose_nash",
    "indef_mod",
    "integreq",
    "liarwhd",
    "morebv",
    "noncvxu2",
    "noncvxun",
    "nondia",
    "nondquar",
    "penalty1",
    "penalty2",
    "penalty3",
    "power",
    "quartc",
    "sbrybnd",
    "schmvett",
    "scosine",
    "sinquad",
    "tointgss",
    "tquartic",
    "tridia",
    "vardim"
]

function get_test_set(n::Int=40)
    """
    Generator of problems listed in public/data/test-problems.json
    """
    return (
        eval(Meta.parse("ADNLPProblems.$pb(n=$n)")) for pb in get_test_set_names()
    )
end

function get_problem(name::String, n::Int=40)
    """
    Get a problem by name
    """
    return eval(Meta.parse("ADNLPProblems.$(name)(n=$n)"))
end

function dump_public_json(filename::String, data::Dict)
    """
    Write json data to "../public/data/filename.json"
    """
    open("../public/$(filename).json", "w") do io
        JSON.print(io, data, 2)
    end
    return nothing
end

function load_public_json_data(filename::String)
    """
    Load json data from "../public/data/filename.json"
    """
    data = Dict()
    open("../public/$(filename).json", "r") do io
        data = JSON.parse(io)
    end
    return data
end
   
function write_to_csv(df::DataFrames.DataFrame, filename::String)
    """
    Write DataFrame to csv file
    """
    CSV.write("../public/$(filename).csv", df)
    return nothing
end

function read_from_csv(filename::String)
    """
    Read DataFrame from csv file
    """
    df = CSV.read("../public/$(filename).csv")
    return df
end

function build_sample_box(problem::ADNLPModel)
    """ 
    Box centered around the minimizer

    The box is built by taking the absolute value of the initial point
    and scaling it by 2. The box is then centered around the initial point.
    The box is used to sample points uniformly around the initial point.
    """
    x0 = problem.meta.x0
    scale_vector = 2 .* abs.(x0)
    scale_vector[scale_vector .<= 1.0] .= 1.0
    return scale_vector
end


function pull_sample(problem::ADNLPModel, box::Vector)
    """
    Pulls a sample uniformly from the box box surrounding x0
    """
    x0_norm = norm(problem.meta.x0)
    if x0_norm < 1.0
        x0_norm = 1.0
    end
    return rand.(Uniform.(-x0_norm, x0_norm), length(problem.meta.x0))
end

function bfgs_linesearch()
    """
    Define the algorithm for the optimization
    """
    return BFGS(;
        alphaguess = LineSearches.InitialStatic(),
        linesearch = LineSearches.HagerZhang(),
        initial_invH = x -> Matrix{eltype(x)}(I, length(x), length(x)),
        manifold = Flat(),
    )
end