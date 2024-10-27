"""
This module automates the extraction, transformation, and loading of optimization programs into various formats.

Globals:
    FILE_PATH: str - The file path of the article this script is based on.
    llm: ChatOpenAI - The OpenAI language model used for generating responses.
"""

import re
from os import getenv
from pypdf import PdfReader
from dotenv import load_dotenv
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.output_parsers import StrOutputParser
import openai as _openai
import logging

load_dotenv()

_openai.api_key = getenv("OPENAI_API_KEY", "")

FILE_PATH = 'references/An-Uncononstrained-Optimization-Test-Functions-Collection.pdf'

llm = ChatOpenAI(model='chatgpt-4o-latest')

def load_pdf() -> str:
    """Loads the PDF file and extracts its text."""
    reader = PdfReader(FILE_PATH)
    return ''.join(page.extract_text() for page in reader.pages)

def clean_text(text: str) -> str:
    """Removes extra spaces, tabs, and empty lines from the text."""
    return re.sub(r'\s+', ' ', text).strip()

def extract_program_list_from_text(text: str) -> str:
    """
    Extracts the list of optimization programs from the text.
    Ignores all text until '2. Unconstrained Optimization Test Functions' appears
    and stops extracting at the 'References' section.
    """
    programs = []
    start = False
    for line in text.split('\n'):
        if start:
            programs.append(line)
        if '2. Unconstrained Optim' in line:
            start = True
            print("Found the start of the optimization programs.")
            continue
        elif 'Refe renc e' in line:
            print("Reached the end of the optimization programs.")
            programs.pop()  # Remove the last line, which is the 'References' line
            break

    return clean_text("\n".join(programs))

from langchain.prompts import SystemMessagePromptTemplate, HumanMessagePromptTemplate

def generate_programs(programs: str) -> list[dict[str, str]]:
    """
    Generates the optimization programs by interacting with the LLM.
    Takes a list of raw program descriptions and returns structured data.
    """
    programs_str = "\n".join(programs)

    # Use 2-tuples for system and user messages as per documentation
    program_latex_generation_template = ChatPromptTemplate.from_messages(
        [
            ("system", "You are a helpful AI assistant that generates optimization programs from raw text."),
            ("human", """Here is a string of unconstrained optimization problems. Please parse the problem name, objective function, and initial iterate: {{{programs_str}}}""")
        ]
    )

    print(f"Program LaTeX generation template created: {program_latex_generation_template}")

    # Use LangChain to invoke the chain and parse the output
    chain = program_latex_generation_template | llm | StrOutputParser()

    # Now invoke the chain with the correct variable
    result = chain.invoke({"programs_str": programs_str})

    print(f"Generated Programs: {result}")

    return result


def main():
    # Load and clean the PDF text
    raw_text = load_pdf()

    # Extract the list of optimization programs
    programs = extract_program_list_from_text(raw_text)

    # Generate the formatted optimization programs
    formatted_programs = generate_programs(programs)

    # Output the formatted programs for debugging
    print(formatted_programs)
