import logging

from ..utils.base_tool import BaseTool
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.smiles import is_smiles


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SmartsCheck(BaseTool):
    __version__ = "0.1.0"
    name = "SmartsCheck"
    func_name = 'check_smarts'
    description = "Check the syntactical validity of a reaction SMART string ([reactant SMILES]>[reagent SMILES]>[product SMILES])."
    categories = ["Molecule", "Reaction"]
    tags = ["SMILES", "SMARTS", "RDKit", "Molecular Information", "Reaction Information"]
    required_envs = []
    text_input_sig = [("smarts", "str", "N/A", "The SMARTS string of a chemical reaction to check.")]
    code_input_sig = [("smarts", "str", "N/A", "The SMARTS string of a chemical reaction to check.")]
    output_sig = [("result", "str", "Description of the validity of the SMARTS string.")]
    examples = [
        {
            'text_input': {'smarts': 'B.C1=CCCCC1.C1=CCCCC1>>B(C1CCCCC1)C1CCCCC1'}, 
            'code_input': {'smarts': 'B.C1=CCCCC1.C1=CCCCC1>>B(C1CCCCC1)C1CCCCC1'}, 
            'output': {'result': 'The reaction SMILES string is valid.'}
        },
    ]

    def _run_base(self, smarts: str) -> str:
        parts = smarts.split('>')
        if len(parts) != 3:
            return f'The SMILES string contains ">", which indicates a reaction, but it is not a valid reaction because it contains {len(parts)} parts instead of 3 (reactants > reagents > products).'
        
        reactants, reagents, products = parts
        reactants = reactants.strip()
        reagents = reagents.strip()
        products = products.strip()

        reactants_valid = is_smiles(reactants)
        if reagents == '':
            reagents_valid = True
        else:
            reagents_valid = is_smiles(reagents)
        products_valid = is_smiles(products)
        
        if reactants_valid and reagents_valid and products_valid:
            return "The reaction SMILES string is valid."
        else:
            invalid_parts = []
            if not reactants_valid:
                invalid_parts.append('reactant(s)')
            if not reagents_valid:
                invalid_parts.append('reagent(s)')
            if not products_valid:
                invalid_parts.append('product(s)')
            return "The reaction SMILES string is invalid. Specifically, the %s is/are not valid SMILES string(s)." % ', '.join(invalid_parts)


if __name__ == "__main__":
    run_mcp_server()
