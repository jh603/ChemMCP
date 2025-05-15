import logging

from ..utils.base_tool import BaseTool, register_mcp_tool
from ..utils.errors import ChemMTKInputError, ChemMTKToolProcessError
from ..utils.smiles import is_smiles
from ..utils.names import smiles2formula
from ..utils.mcp_app import mcp_instance, run_mcp_server


logger = logging.getLogger(__name__)


@register_mcp_tool(mcp_instance)
class Smiles2Formula(BaseTool):
    __version__ = "0.1.0"
    name = "SMILES2Formula"
    func_name = 'convert_smiles_to_formula'
    description = "Convert SMILES to molecular formula."
    categories = ["Molecule"]
    tags = ["SMILES", "RDKit"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    output_sig = [('formula', 'str', 'Molecular formula of the molecule')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'formula': 'C2H6O'}},
    ]

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMTKInputError("The input is not a valid SMILES string.")
        
        try:
            formula = smiles2formula(smiles)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise ChemMTKToolProcessError("Failed to process the SMILES string.") from e
        
        return formula


if __name__ == "__main__":
    run_mcp_server()
