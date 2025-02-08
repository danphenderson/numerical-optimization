module Common

using LinearAlgebra
using OptimizationProblems, 
using JSON
using DataFrames
using Plots
using CSV

export get_test_set, dump_public_json, load_public_json_data, write_to_csv, read_from_csv, get_problem, get_test_set_names

get_test_set_names() = load_public_json_data("test-problems")

function get_test_set(n::Int=40)
    """
    Generator of problems listed in public/data/test-problems.json
    """
    return (
        eval(Meta.parse("ADNLPProblems.$pb(n=$n)")) for pb in load_public_json_data("test-problems")
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

end

