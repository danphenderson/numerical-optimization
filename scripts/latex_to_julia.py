"""
Script to automate the transformation of LaTeX
objective functions to Julia functions.
"""

import re
from os import getenv
from dotenv import load_dotenv
from pathlib import Path
from langsmith import expect
import openai as _openai
import logging
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

load_dotenv()

_openai.api_key = getenv("OPENAI_API_KEY", "")

_PATH = Path(__file__).resolve()

_TEST_SET_PATH = _PATH.parent.parent / Path("test-set") / Path("programs.tex")

_OUTPUT_PATH = _PATH.parent.parent / Path("output") / Path("julia_functions.jl")

llm = ChatOpenAI(model='chatgpt-4o-latest')

# Initialize the OpenAI API
def initialize_openai():
    if not _openai.api_key:
        raise ValueError("OpenAI API key not found. Please set the OPENAI_API_KEY environment variable.")
    return _openai

def load_latex_enumeration(file_path: Path) -> str:
    """Load LaTeX content from the specified file path."""
    try:
        with open(file_path, 'r') as file:
            return file.read()
    except FileNotFoundError:
        logging.error(f"File not found: {file_path}")
    except Exception as e:
        logging.error(f"Error loading LaTeX content: {e}")
    return ""

def parse_latex(latex_string) -> list[str]:
    """Parses each item in a latex enumeration environment."""
    # Extract the content within \begin{enumerate} and \end{enumerate}
    enumerate_pattern = r'\\begin{enumerate}(.*?)\\end{enumerate}'
    enumerate_content = re.findall(enumerate_pattern, latex_string, re.DOTALL)

    if not enumerate_content:
        logging.error("No 'enumerate' environment found in the LaTeX content.")
        return []

    # Extract each \item content
    item_pattern = r'\\item\s*(.*?)\s*(?=\\item|\\end{enumerate})'
    items = re.findall(item_pattern, enumerate_content[0], re.DOTALL)

    parsed_items = [item.strip() for item in items if item.strip()]
    return parsed_items

def convert_parsed_latex_to_julia(parsed_latex) -> list[str]:
    """Returns a list of Julia functions as strings.

    Input `parsed_latex` is a list of strings, where
    each string represents an objective function of an
    unconstrained optimization program.
    """


    latex_to_julia_template = ChatPromptTemplate.from_messages(
        [
            (
                "system",
                """You are a helpful AI assistant that converts LaTeX optimization programs to Julia functions. Specifically, you are given the LaTeX objective function, initial iterate, and problem name for an unconstrained optimization program.""",
            ),
            (
                "human",
                """For example, consider this LaTeX program:

```tex
\\textbf{{Generalized Rosenbrock function:}}
    \\[f(x) = \\sum_{{i=1}}^{{n-1}} \\left( c(x_{{i+1}} - x_i^2)^2 + (1 - x_i)^2 \\right),\\]
    \\[x_0 = [-1.2, 1, -1.2, 1, \\dots, -1.2, 1], \\quad c = 100.\\]
```

The expected output should look like this:

```julia
f = x -> begin
    N = lastindex(x)
    return 100sum((x[i+1] - x[i]^2)^2 for i = 1:N-1) + sum((x[i] - 1.0)^2 for i = 1:N-1)
end

init = (n::Int=500) -> begin
    x0 = [j/(n+1) for j in 1:n]
    return n, x0, "Generalized Rosenbrock"
end
```

You should only return the ```julia``` code block when prompted with a LateX optimization program. Return nothing more than the Julia code block.
"""),
            ("ai", "Great! Can you provide me with one more example?"),
            ("human", """Another example, consider this LaTeX program:

```tex
'\\textbf{{Extended Penalty function:}}
    \\[f(x) =  \\sum_{{i=1}}^{{n-1}} (x_i - 1)^2 + \\left(\\sum_{{j=1}}^{{n}}  x_j^2 - 0.25 \\right)^2,\\]
    \\[x_0 = [1, 2, \\dots, n].\\]'

The expected output should look like this:

```julia
f = x -> begin
    N = lastindex(x)
    return sum((x[i] - 1.0)^2 for i = 1:N) + (sum(x[j]^2 for j = 1:N) - 0.25)^2
end

init = (n::Int=500) -> begin
    x0 = [j for j in 1:n]
    return n, x0, "Extended Penalty"
end
```"""),
            ("ai", "Great! provide me with with LaTeX programs and I will generate the corresponding Julia code blocks."),
            ("human", "{item}"),
        ]
    )

    julia_functions = []

    chain = latex_to_julia_template | llm | StrOutputParser()
    for idx, item in enumerate(parsed_latex[2:]):
        print(f"Processing LaTeX item {idx + 1}...")
        output = chain.invoke({"item": item})
        julia_functions.append(output)

    return julia_functions


def main():
    raw_latex = load_latex_enumeration(_TEST_SET_PATH)

    if not raw_latex:
        print("Failed to load LaTeX content.")
        return

    parsed_latex = parse_latex(raw_latex)

    if not parsed_latex:
        print("No objective functions found in the LaTeX content.")
        return

    julia_functions = convert_parsed_latex_to_julia(parsed_latex)

    if not julia_functions:
        print("No Julia functions were generated.")
        return

    with open(_OUTPUT_PATH, 'w') as file:
        for idx, julia_function in enumerate(julia_functions):
            file.write(julia_function)
            file.write("\n\n")

    print(f"Julia functions written to {_OUTPUT_PATH}")

if __name__ == "__main__":
    main()
