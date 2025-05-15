import logging

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMTKSearchFailError
from ..tool_utils.names import pubchem_name2smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class Name2Smiles(BaseTool):
    __version__ = "0.1.0"
    name = "Name2Smiles"
    func_name = 'convert_chemical_name_to_smiles'
    description = "Convert chemical name to SMILES string."
    categories = ["Molecule"]
    tags = ["SMILES", "RDKit", "PubChem"]
    required_envs = []
    code_input_sig = [('name', 'str', 'N/A', 'Chemical name of the molecule')]
    text_input_sig = [('name', 'str', 'N/A', 'Chemical name of the molecule')]
    output_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    examples = [
        {'code_input': {'name': 'aspirin'}, 'text_input': {'name': 'aspirin'}, 'output': {'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'}},
    ]

    def _run_base(self, name: str) -> str:
        try:
            smi = pubchem_name2smiles(name)
            logger.debug("Looking up PubChem succeeded.")
        except ChemMTKSearchFailError as e:
            logger.debug("Looking up PubChem failed.")
            raise e
        
        return smi


if __name__ == "__main__":
    run_mcp_server()
