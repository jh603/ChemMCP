import argparse
import requests
from copy import deepcopy
from dataclasses import dataclass
from typing import Optional, Literal

import pubchempy as pcp

from .utils.base_tool import BaseTool, register_mcp_tool
from .utils.errors import *
from .utils.smiles import is_smiles
from .utils.pubchem import pubchem_iupac2cid, pubchem_name2cid
from .utils.llm import llm_completion
from .utils.mcp_app import mcp_instance


QA_SYSTEM_PROMPT = "You are an expert chemist. You will be given the PubChem page about a molecule/compound, and your task is to answer the question based on the information of the page. Your answer should be accurate and concise, and contain all the information necessary to answer the question."


unuseful_section_names = {
    "Structures": None,
    "Chemical Safety": None,
    "Names and Identifiers": {
        "Other Identifiers": None,
        "Synonyms": None,
        "Create Date": None,
        "Modify Date": None,
    },
    "Chemical and Physical Properties": {
        "SpringerMaterials Properties": None,
    },
    "Spectral Information": None,
    "Related Records": None,
    "Chemical Vendors": None,
    "Drug and Medication Information": {
        "WHO Essential Medicines": None,
        "FDA Approved Drugs": None,
        "FDA Orange Book": None,
        "FDA National Drug Code Directory": None,
        "FDA Green Book": None,
        "Drug Labels": None,
        "Clinical Trials": None,
        # "Therapeutic Uses": None,
        # "Drug Warnings": None,
        # "Drug Idiosyncrasies": None,
        # "Reported Fatal Dose": None,
        # "Maximum Drug Dose": None,
    },
    "Pharmacology and Biochemistry": None,
    "Use and Manufacturing": None,
    "Identification": None,
    "Literature": None,
    "Patents": None,
    "Interactions and Pathways": None,
    "Biological Test Results": None,
    "Classification": None,
    "Taxonomy": None,
}


@dataclass
class Information:
    information_item: dict

    @classmethod
    def construct(cls, data):
        information = cls(data)
        return information

    def generate_text(self, display_controls=None):
        data = self.information_item
        value = data['Value']
        text = ""
        if 'StringWithMarkup' in value:
            strings = value['StringWithMarkup']
            for item in strings:
                tmp_text = item['String']
                tmp_unit = (" " + item['Unit']) if 'Unit' in item else ''
                text += tmp_text + tmp_unit + '\n'
        elif 'Number' in value:
            if 'Name' in value:
                name = value['Name']
                text += name + ": "
            strings = value['Number']
            strings = [str(item) for item in strings]
            text += ', '.join(strings)
            tmp_unit = (" " + value['Unit']) if 'Unit' in value else ''
            text += tmp_unit + '\n'
        
        if text.strip() == "":
            return None
        text = text
        return text


@dataclass
class Section:
    level: int
    title: str
    description: str
    display_controls: dict
    information_list: Optional[list] = None
    subsection_list: Optional[list] = None

    @classmethod
    def construct(cls, data, level=1):
        title = data['TOCHeading']
        description = data['Description'] if 'Description' in data else None
        display_controls = data['DisplayControls'] if 'DisplayControls' in data else None

        section = cls(level=level, title=title, description=description, display_controls=display_controls)

        if 'Information' in data:
            information_list = []
            for information_data in data['Information']:
                information = Information.construct(information_data)
                information_list.append(information)
            section.information_list = information_list
        
        if 'Section' in data:
            subsection_list = []
            for subsection_data in data['Section']:
                subsection = Section.construct(subsection_data, level=level + 1)
                subsection_list.append(subsection)
            section.subsection_list = subsection_list

        return section
        
    def generate_text(self, indices=None) -> str:
        # if self.display_controls is not None and "HideThisSection" in self.display_controls and self.display_controls["HideThisSection"] is True:
        #     return None
        
        title_text = '#' * self.level + ' ' + (('.'.join(indices) + ' ') if len(indices) > 0 else '') + self.title + '\n'
        if self.description is not None:
            title_text += 'Section Description: ' + self.description
        title_text += '\n\n'

        if indices is None:
            indices = tuple()

        content_text = ""

        if self.information_list is not None:
            for information in self.information_list:
                tmp_text = information.generate_text(self.display_controls)
                if tmp_text is not None:
                    content_text += tmp_text

        if self.subsection_list is not None:
            idx = 1
            for subsection in self.subsection_list:
                tmp_text = subsection.generate_text(indices + (str(idx),))
                if tmp_text is not None:
                    idx += 1
                    content_text += tmp_text
        
        if content_text.strip() == "":
            return None

        text = title_text + content_text + '\n\n'

        return text


@dataclass
class PubchemStructuredDoc:
    doc_data = []

    @classmethod
    def construct(cls, sections):
        doc = PubchemStructuredDoc()
        section_list = []
        
        for section_data in sections:
            section = Section.construct(section_data)
            section_list.append(section)
        
        doc.doc_data = section_list

        return doc

    def generate_text(self) -> str:
        text = ""
        idx = 1
        for section in self.doc_data:
            tmp_text = section.generate_text(indices=(str(idx),))
            if tmp_text is not None:
                text += tmp_text
                idx += 1
        return text


class PubchemSearch(BaseTool):
    name = "PubchemSearch"
    func_name = 'search_pubchem'
    description = "Search for molecule/compound information on PubChem, one of the most comprehensive database of chemical molecules and their activities."
    text_input_sig = [("representation_name_and_representation", "str", "The representation name and representation of the molecule/compound, e.g., \"SMILES: <SMILES>\", \"IUPAC: <IUPAC name>\", or \"Name: <common name>\".")]
    code_input_sig = [("representation_name", "str", "The representation name, can be \"SMILES\", \"IUPAC\", or \"Name\" (chemical's common name)."), ("representation", "str", "The representation of the molecule/compound, corresponding to the representation_name used.")]
    output_sig = [("compound_doc", "str", "The document of the molecule/compound in a markdown format.")]
    examples = [
        {'code_input': {'representation_name': 'SMILES', 'representation': 'CCO'}, 'text_input': {'representation_name_and_representation': 'SMILES: CCO'}, 'output': {'compound_doc': '# 1 Names and Identifiers\nSection Description: Chemical names, synonyms, identifiers, and descriptors.\n\n## 1.1 Record Description\nSection Description: Summary Information\n\nEthanol with a small amount of an adulterant added so as to be unfit for use as a beverage. [...]'}},
        
        {'code_input': {'representation_name': 'IUPAC', 'representation': 'ethanol'}, 'text_input': {'representation_name_and_representation': 'IUPAC: ethanol'}, 'output': {'compound_doc': '# 1 Names and Identifiers\nSection Description: Chemical names, synonyms, identifiers, and descriptors.\n\n## 1.1 Record Description\nSection Description: Summary Information\n\nEthanol with a small amount of an adulterant added so as to be unfit for use as a beverage. [...]'}},

        {'code_input': {'representation_name': 'Name', 'representation': 'alcohol'}, 'text_input': {'representation_name_and_representation': 'Name: alcohol'}, 'output': {'compound_doc': '# 1 Names and Identifiers\nSection Description: Chemical names, synonyms, identifiers, and descriptors.\n\n## 1.1 Record Description\nSection Description: Summary Information\n\nEthanol with a small amount of an adulterant added so as to be unfit for use as a beverage. [...]'}},
    ]
    
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/JSON/'

    def _run_text(self, representation_name_and_representation: str) -> str:
        try:
            namespace, identifier = representation_name_and_representation.split(':')
            namespace = namespace.strip()
            identifier = identifier.strip()
            if namespace.lower() in ('smiles',):
                namespace = 'smiles'
            elif namespace.lower() in ('iupac', 'iupac name'):
                namespace = 'iupac'
            elif namespace.lower() in ('name', 'common name'):
                namespace = 'name'
            elif namespace == '':
                raise ChemMTKInputError('Empty representation name.')
            else:
                raise ChemMTKInputError('The representation name \"%s\" is not supported. Please use \"SMILES\", \"IUPAC\", or \"Name\".' % namespace)
        except (ChemMTKInputError, ValueError) as e:
            raise ChemMTKInputError("The input is not in a correct format: %s If searching with SMILES, please input \"SMILES: <SMILES of the molecule/compound>\"; if searching with IUPAC name, please input \"IUPAC: <IUPAC name of the molecule/compound>\"; if searching with common name, please input \"Name: <common name of the molecule/compound>\"." % str(e))
        r = self._run_base(namespace, identifier)
        return r
    
    def _run_base(self, representation_name: Literal["smiles", 'iupac', 'name'], representation: str) -> str:
        cid = self._search_cid(representation_name, representation)
        return self.get_cid_doc_text(cid)
        
    def get_cid_doc_text(self, cid):
        data = self.get_data(cid)

        try:
            sections = self.remove_unuseful_sections(data['Record']['Section'])
        except KeyError:
            print(data)
            print('cid: ', cid)
            raise

        doc = self.construct_doc_text(sections)
        
        return doc

    def _search_cid(self, namespace, identifier):
        if namespace == 'smiles' and not is_smiles(identifier):
            raise ChemMTKInputError('The input SMILES is invalid. Please double-check. Note that you should input only one molecule/compound at a time.')
        
        if namespace == 'iupac':
            cid = pubchem_iupac2cid(identifier)
        elif namespace == 'smiles':
            try:
                c = pcp.get_compounds(identifier, namespace=namespace)
            except pcp.BadRequestError:
                raise ChemMTKSearchFailError("Error occurred while searching for the molecule/compound on PubChem. Please try other tools or double check your input.")
            if len(c) >= 1:
                c = c[0]
            else:
                raise ChemMTKSearchFailError("Could not find a matched molecule/compound on PubChem. Please double check your input and search for one molecule/compound at a time, or use its another identifier (e.g., IUPAC name or common name) for the search.")
            cid = c.cid
        else:
            cid = pubchem_name2cid(identifier)
        
        if cid is None:
            raise ChemMTKSearchFailError("Could not find a matched molecule/compound on PubChem. Please double check your input and search for one molecule/compound at a time, or use its another identifier for the search.")

        return cid
    
    @staticmethod
    def get_data(cid):
        url = PubchemSearch.url.format(cid)
        data = requests.get(url).json()
        return data
    
    @staticmethod
    def construct_doc_text(sections):
        doc = PubchemStructuredDoc.construct(sections)
        text = doc.generate_text()
        return text
    
    @staticmethod
    def remove_unuseful_sections(sections):
        sections = deepcopy(sections)

        new_sections = []
        for section in sections:

            section_title = section['TOCHeading']
            if section_title in unuseful_section_names and unuseful_section_names[section_title] is None:
                continue

            if 'Section' in section:
                subsection_list = section['Section']
                new_subsection_list = []
                for subsection in subsection_list:
                    subsection_title = subsection['TOCHeading']
                    if section_title in unuseful_section_names and subsection_title in unuseful_section_names[section_title] and unuseful_section_names[section_title][subsection_title] is None:
                        continue
                    new_subsection_list.append(subsection)
                
                if len(new_subsection_list) == 0:
                    continue

                section['Section'] = new_subsection_list

            new_sections.append(section)
        
        return new_sections
    

@register_mcp_tool(mcp_instance)
class PubchemSearchQA(BaseTool):
    name = "PubchemSearchQA"
    func_name = 'search_pubchem_qa'
    description = "Answer questions about molecules/compounds based on the information from PubChem, one of the most comprehensive database of chemical molecules and their activities. Input \"representation name: representation\" (e.g., \"SMILES: <SMILES>\", \"IUPAC: <IUPAC name>\", or \"Name: <common name>\", one at a time), followed by \"Question: <your question about the molecule/compound>\", returns the related information."
    text_input_sig = [("representation_name_and_representation", "str", "The representation name and representation of the molecule/compound, e.g., \"SMILES: <SMILES>\", \"IUPAC: <IUPAC name>\", or \"Name: <common name>\". Followed by \"Question: <your question about the molecule/compound>\".")]
    code_input_sig = [("representation_name", "str", "The representation name, can be \"smiles\", \"iupac\", or \"name\" (chemical's common name)."), ("representation", "str", "The representation of the molecule/compound, corresponding to the representation_name used."), ("question", "str", "The question about the molecule/compound.")]
    output_sig = [("answer", "str", "The answer to the question based on the PubChem page.")]
    examples = [
        {'code_input': {'representation_name': 'SMILES', 'representation': 'CCO', "question": "What properties do this molecule have?"}, 'text_input': {'representation_name_and_representation_and_question': 'SMILES: CCO   Questions: What properties do this molecule have?'}, 'output': {'answer': 'This molecule has the following properties: [...]'}},
        
        {'code_input': {'representation_name': 'IUPAC', 'representation': 'ethanol', "question": "What properties do this molecule have?"}, 'text_input': {'representation_name_and_representation_and_question': 'IUPAC: ethanol   Questions: What properties do this molecule have?'}, 'output': {'answer': 'This molecule has the following properties: [...]'}},

        {'code_input': {'representation_name': 'Name', 'representation': 'alcohol', "question": "What properties do this molecule have?"}, 'text_input': {'representation_name_and_representation_and_question': 'Name: alcohol   Questions: What properties do this molecule have?'}, 'output': {'answer': 'This molecule has the following properties: [...]'}},
    ]

    def __init__(self, llm_model=None, init=True, interface='text') -> None:
        super().__init__(init, interface)
        self.llm_model = llm_model
        self.pubchem_search = PubchemSearch(init=init, interface='code')

    def _run_text(self, representation_name_and_representation_and_question: str) -> str:
        if 'Question:' not in representation_name_and_representation_and_question:
            raise ChemMTKInputError("The input is not in a correct format. Please input the molecule/compound representation followed by the question about the molecule/compound. An example: \"SMILES: <SMILES of the molecule/compound> Question: <your question about the molecule/compound>\".")
        representation_name_and_representation_and_question, question = representation_name_and_representation_and_question.split('Question:')
        representation_name_and_representation_and_question = representation_name_and_representation_and_question.strip()
        question = question.strip()

        try:
            namespace, identifier = representation_name_and_representation_and_question.split(':')
            namespace = namespace.strip()
            identifier = identifier.strip()
            if namespace.lower() in ('smiles',):
                namespace = 'smiles'
            elif namespace.lower() in ('iupac', 'iupac name'):
                namespace = 'iupac'
            elif namespace.lower() in ('name', 'common name'):
                namespace = 'name'
            elif namespace == '':
                raise ChemMTKInputError('Empty representation name.')
            else:
                raise ChemMTKInputError('The representation name \"%s\" is not supported. Please use \"SMILES\", \"IUPAC\", or \"Name\".' % namespace)
        except (ChemMTKInputError, ValueError) as e:
            raise ChemMTKInputError("The input is not in a correct format: %s If searching with SMILES, please input \"SMILES: <SMILES of the molecule/compound>\"; if searching with IUPAC name, please input \"IUPAC: <IUPAC name of the molecule/compound>\"; if searching with common name, please input \"Name: <common name of the molecule/compound>\". After that, append your question about the molecule/compound as \"Question: <your question>\"." % str(e))
        r = self._run_base(namespace, identifier, question)
        return r
    
    def _run_base(self, representation_name: Literal["SMILES", "IUPAC", "Name"], representation: str, question: str):
        doc = self.pubchem_search.run_code(representation_name, representation)
        conversation = [
            {'role': 'system', 'content': QA_SYSTEM_PROMPT},
            {'role': 'user', 'content': doc + '\n\n\n\nQuestion: ' + question},
        ]
        r = llm_completion(messages=conversation, llm_model=self.llm_model)
        r = r.choices[0].message.content.strip()
        return r


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the MCP server.")
    parser.add_argument('--sse', action='store_true', help="Run the server with SSE (Server-Sent Events) support.")
    args = parser.parse_args()

    if args.sse:
        # build a Starlette/uvicorn app
        app = mcp_instance.sse_app()
        import uvicorn
        uvicorn.run(app, host="127.0.0.1", port=8001)
    else:
        # Run the MCP server with standard input/output
        mcp_instance.run(transport='stdio')
