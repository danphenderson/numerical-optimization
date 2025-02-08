"""
Script to automate the transformation of LaTeX
objective functions to Julia functions.
"""

import re
from os import getenv
from dotenv import load_dotenv
from pathlib import Path

import openai as _openai
import logging

from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser

load_dotenv()

_openai.api_key = getenv("OPENAI_API_KEY", "")

_PATH = Path(__file__).resolve()

_TEST_SET_PATH = _PATH.parent.parent / Path("test-set") / Path("filtered") / Path("interesting-programs.tex")

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
            ("system", "Convert the following LaTeX program to Julia code:"),
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
