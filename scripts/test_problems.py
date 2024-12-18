"""
Module converts a subset of `OptimizationProblems.jl` into mathematica.

We have 80 problems filtered from OptimizationProblems.jl (as shown in
documentation). We start by identifying the location within the OptimizationProblems.jl
repository where the problems are stored. Then, we parse each problems source code
declaration and hand it to the OpenAI API to convert the Julia code into Mathematica code.

Running the script will output the Mathematica code for each problem in the console and store
the output in a file called `test-problems-out.txt`.

TODO:
- It appears that we are only finding 71 of the 80 models.
- All of our 80 models are variable, but some have constraints on the variables.
  This is not currently handled in the mathematica code, as it is in their Julia source
  code declaration. We should either add constraints to the Mathematica code or remove
  such problems from the test set.
"""

from pathlib import Path
from dotenv import load_dotenv
from os import getenv
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser
import openai as _openai
import logging

load_dotenv()

_openai.api_key = getenv("OPENAI_API_KEY", "")

llm = ChatOpenAI(model='chatgpt-4o-latest')

MODELS = [
    # "cosine",
    # "arglinb","arglinc","argtrig","arwhead","bdqrtic","bearing",
    # "brownal","broyden3d","broydn7d","chainwoo","chnrosnb_mod","clplatea",
    # "clplateb","clplatec","cosine","cragglvy","cragglvy2","curly","curly10","curly20",
    # "curly30", "dixon3dq","dqdrtic","dqrtic",
    # "edensch","engval1","errinros_mod","extrosnb","fletcbv2","fletcbv3_mod","fletchcr",
    # "genhumps","genrose_nash","indef_mod","integreq","liarwhd","morebv",
    # "ncb20","ncb20b","noncvxu2","noncvxun","nondia","nondquar","penalty1","penalty2","penalty3",
    # "power","quartc","sbrybnd","schmvett","scosine","sinquad","sparsine","sparsqur",
    # "srosenbr","tointgss","tquartic","tridia","vardim","woods"

    # EXCLUDED PROBLEMS:
    # The following problems are given as examples:
    # "genrose", "arglina", "eg2", "freuroth", "brybnd",    # Converted Manually
    # The following problems have constraints on the dimension:
    # "dixmaane", "dixmaanf","dixmaang","dixmaanh","dixmaani","dixmaanj",   # Multiple of 3
    # "dixmaank","dixmaanl","dixmaanm","dixmaann","dixmaano","dixmaanp",    # Multiple of 3
    # "nzf1",       # Must be greater than 26
    # "powellsg"    # Multiple of 4
    # "spmsrtls"    # Ref: https://github.com/JuliaSmoothOptimizers/OptimizationProblems.jl/issues/354
]

FILES = [file for file in (Path(".") / "ADNLPProblems").iterdir()]

OUTPUT_PATH = Path(".") / "test-problems-out.txt"


def clean_normalized_path(path: Path) -> str:
    """Used to match the file corresponding to all of our models."""
    return path.name.lower().replace("_", "").replace("-", "").lower()

def models_files_dict()-> dict:
    """
    For each model in models, we search for the corresponding
    file in the ADNLPModels directory.

    Note, the files are stored in the files list and the
    models are stored in the models list.

    FIXME
    - Only finding 71 of the 80 models.
    """
    model_to_file = {}
    for model in MODELS:
        for file in FILES:
            if model in clean_normalized_path(file):
                model_to_file[model] = file
            if model in file.name.lower():
                model_to_file[model] = file
    return model_to_file


def load_julia_code(file_path: Path) -> str:
    """
    Load the Julia code from the specified file path.
    """
    try:
        with open(file_path, 'r') as file:
            return file.read()
    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
    except Exception as e:
        logging.error(f"Error loading Julia code: {e}")
    return ""

def julia_to_mathematica():
    """
    Converts the Julia code to LaTeX equations.
    """
    model_to_file = models_files_dict()

    prompt = ChatPromptTemplate.from_messages(
        [
            (
                "system",
                "You are an expert extraction algorithm. "
                "Your task is to extract the name, objective function definition, and x0 associated "
                "with a Julia ADNLPModel instance and return the corresponding Mathematica code. "
                "Specifically, return ONLY the corresponding association declared using the following helper function. "
                """`CreateTestFunction[name_, def_, x0_] := <|"Name" -> name, "Definition" -> def, "x0" -> x0|>;`"""
            ),
            (
                "ai",
                "Great! Please provide an example Julia code input and the expected Mathematica code output."
            ),
            (
                "human",
                """Example Julia code input:
export genrose, rosenbrock

function genrose(; n::Int = default_nvar, type::Type{{T}} = Float64, kwargs...) where {{T}}
  n < 2 && @warn("genrose: number of variables must be ≥ 2")
  n = max(2, n)
  function f(x; n = length(x))
    return 1 +
           100 * sum((x[i + 1] - x[i]^2)^2 for i = 1:(n - 1)) +
           sum((x[i] - 1)^2 for i = 1:(n - 1))
  end
  x0 = T.([i / (n + 1) for i = 1:n])
  return ADNLPModels.ADNLPModel(f, x0, name = "genrose"; kwargs...)
end

rosenbrock(args...; kwargs...) = genrose(args..., n = 2; delete!(Dict(kwargs), :n)...)


Expected Mathematica code output:
CreateTestFunction[
    "genrose",
    Function[x, Module[{{n = Length[x]}},
        Sum[100*(x[[i + 1]] - x[[i]]^2)^2 + (1.0 - x[[i]])^2, {{i, n - 1}}]
    ]],
    Function[n, Table[If[OddQ[i], -1.2, 1], {{i, n}}]]
]
"""
            ),
            ("ai", "Can you provide me with another example Julia code input and the expected Mathematica code output?"),
            ("human", """Example Julia code input:
export arglina

function arglina(; use_nls::Bool = false, kwargs...)
  model = use_nls ? :nls : :nlp
  return arglina(Val(model); kwargs...)
end

function arglina(::Val{{:nlp}}; n::Int = default_nvar, type::Type{{T}} = Float64, kwargs...) where {{T}}
  function f(x; n = length(x))
    m = 2 * n
    sj = sum(x[j] for j = 1:n)
    return 1 // 2 * sum((x[i] - 2 // m * sj - 1)^2 for i = 1:n) +
           1 // 2 * sum((-2 // m * sj - 1)^2 for i = (n + 1):m)
  end
  x0 = ones(T, n)
  return ADNLPModels.ADNLPModel(f, x0, name = "arglina"; kwargs...)
end

function arglina(::Val{{:nls}}; n::Int = default_nvar, type::Type{{T}} = Float64, kwargs...) where {{T}}
  function F!(r, x)
    m = 2 * n
    sj = sum(x[j] for j = 1:n)
    for i = 1:n
      r[i] = x[i] - 2 // m * sj - 1
      r[i + n] = -2 // m * sj - 1
    end
    return r
  end
  x0 = ones(T, n)
  return ADNLPModels.ADNLSModel!(F!, x0, 2 * n, name = "arglina-nls"; kwargs...)
end


Expected Mathematica code output:
CreateTestFunction[
    "arglina",
    Function[x, Module[{{n = Length[x], m, sj}},
        m = 2*n;
        sj = Total[x];
        (1/2)*Sum[(x[[i]] - (2/m)*sj - 1)^2, {{i, n}}] +
        (1/2)*Sum[(-(2/m)*sj - 1)^2, {{i, n+1, m}}]
    ]],
    Function[n, ConstantArray[1, n]]
]"""
            ),
            ("ai", "Can you provide me with another example Julia code input and the expected Mathematica code output?"),
            ("human", """Example Julia code input:
export eg2

function eg2(; n::Int = default_nvar, type::Type{{T}} = Float64, kwargs...) where {{T}}
  n < 2 && @warn("eg2: number of variables must be ≥ 2")
  n = max(2, n)
  function f(x; n = length(x))
    sum(sin(x[1] + x[i]^2 - 1) for i = 1:(n - 1)) + sin(x[n]^2) / 2
  end
  x0 = zeros(T, n)
  return ADNLPModels.ADNLPModel(f, x0, name = "eg2"; kwargs...)
end


Expected Mathematica code output:
CreateTestFunction[
    "eg2",
    Function[x,
    Module[{{n = Length[x]}},
        Sum[Sin[x[[1]] + x[[i]]^2 - 1], {{i, 1, n - 1}}] + (1/2)*Sin[x[[n]]^2]
    ]],
    Function[n, ConstantArray[0, n]]
]"""
            ),
            ("ai", "Can you provide me with another example Julia code input and the expected Mathematica code output?"),
            ("human", """Example Julia code input:
export freuroth

function freuroth(; use_nls::Bool = false, kwargs...)
  model = use_nls ? :nls : :nlp
  return freuroth(Val(model); kwargs...)
end

function freuroth(::Val{{:nlp}}; n::Int = default_nvar, type::Type{{T}} = Float64, kwargs...) where {{T}}
  n < 2 && @warn("freuroth: number of variables must be ≥ 2")
  n = max(2, n)
  function f(x; n = length(x))
    return 1 // 2 *
           sum(((5 - x[i + 1]) * x[i + 1]^2 + x[i] - 2 * x[i + 1] - 13)^2 for i = 1:(n - 1)) +
           1 // 2 *
           sum(((1 + x[i + 1]) * x[i + 1]^2 + x[i] - 14 * x[i + 1] - 29)^2 for i = 1:(n - 1))
  end

  x0 = zeros(T, n)
  x0[1] = one(T) / 2
  x0[2] = -2 * one(T)
  return ADNLPModels.ADNLPModel(f, x0, name = "freuroth"; kwargs...)
end

function freuroth(::Val{{:nls}}; n::Int = default_nvar, type::Type{{T}} = Float64, kwargs...) where {{T}}
  n < 2 && @warn("freuroth: number of variables must be ≥ 2")
  n = max(2, n)
  function F!(r, x; n = length(x))
    for i = 1:(n - 1)
      r[i] = (5 - x[i + 1]) * x[i + 1]^2 + x[i] - 2 * x[i + 1] - 13
      r[i + n - 1] = (1 + x[i + 1]) * x[i + 1]^2 + x[i] - 14 * x[i + 1] - 29
    end
    return r
  end

  x0 = zeros(T, n)
  x0[1] = one(T) / 2
  x0[2] = -2 * one(T)
  return ADNLPModels.ADNLSModel!(F!, x0, 2 * (n - 1), name = "freuroth-nls"; kwargs...)
end


Expected Mathematica code output:
CreateTestFunction[
    "freuroth",
    Function[x,
    Module[{{n = Length[x]}},
        (1/2) * Sum[((5 - x[[i + 1]])*x[[i + 1]]^2 + x[[i]] - 2*x[[i + 1]] - 13)^2, {{i, 1, n - 1}}] +
        (1/2) * Sum[((1 + x[[i + 1]])*x[[i + 1]]^2 + x[[i]] - 14*x[[i + 1]] - 29)^2, {{i, 1, n - 1}}]
    ]],
    Function[n,
    Module[{{x0 = ConstantArray[0, n]}},
        x0[[1]] = 1/2;
        x0[[2]] = -2;
        x0
    ]]
]"""
            ),
            ("ai", "Can you provide me with another example Julia code input and the expected Mathematica code output?"),
            ("human", """Example Julia code input:
export brybnd

function brybnd(; use_nls::Bool = false, kwargs...)
  model = use_nls ? :nls : :nlp
  return brybnd(Val(model); kwargs...)
end

function brybnd(::Val{{:nlp}}; n::Int = default_nvar, type::Type{{T}} = Float64, kwargs...) where {{T}}
  function f(x; n = length(x))
    ml = 5
    mu = 1
    return 1 // 2 * sum(
      (
        x[i] * (2 + 5 * x[i]^2) + 1 -
        sum(x[j] * (1 + x[j]) for j = max(1, i - ml):min(n, i + mu) if j != i)
      )^2 for i = 1:n
    )
  end
  x0 = -ones(T, n)
  return ADNLPModels.ADNLPModel(f, x0, name = "brybnd"; kwargs...)
end

function brybnd(::Val{{:nls}}; n::Int = default_nvar, type::Type{{T}} = Float64, kwargs...) where {{T}}
  function F!(r, x; n = length(x))
    ml = 5
    mu = 1
    for i = 1:n
      r[i] =
        x[i] * (2 + 5 * x[i]^2) + 1 -
        sum(x[j] * (1 + x[j]) for j = max(1, i - ml):min(n, i + mu) if j != i)
    end
    return r
  end
  x0 = -ones(T, n)
  return ADNLPModels.ADNLSModel!(F!, x0, n, name = "brybnd-nls"; kwargs...)
end


Expected Mathematica code output:
CreateTestFunction[
    "brybnd",
    Function[x,
    Module[{{n = Length[x], ml = 5, mu = 1}},
        (1/2) * Sum[
            (
                x[[i]] * (2 + 5 * x[[i]]^2) + 1 -
                Sum[If[j != i, x[[j]] * (1 + x[[j]]), 0], {{j, Max[1, i - ml], Min[n, i + mu]}}]
            )^2,
            {{i, 1, n}}
        ]
    ]],
    Function[n, -ConstantArray[1, n]]
]"""
            ),
            (
                "ai",
                "Great! I am ready to begin extraction. "
                "When prompted with the next Julia code input I will perform "
                "my extraction algorithm and return the corresponding Mathematica code as output."),
            (
                "human",
                "Your response should contain ONLY the corresponding Mathematica code. "
                "Here is the next Julia code input: {item}"
            ),
        ]
    )

    res = {}

    for model, file_path in model_to_file.items():
        julia_code = load_julia_code(file_path)

        if not julia_code:
            continue

        chain = prompt | llm | StrOutputParser()
        try:
            output = chain.invoke({"item": julia_code})
            # If output is empty or None, skip
            if not output or not output.strip():
                continue
            print(output)
            res[model] = output
        except Exception as e:
            logging.error(f"Extraction failed for {model}: {e}")
            # Skip this model on failure
            continue
    return res
