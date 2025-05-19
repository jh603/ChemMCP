import os
from typing import Optional

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPApiNotFoundError, ChemMCPInputError
from ..tool_utils.chemspace import ChemSpace
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server


@ChemMCPManager.register_tool
class MoleculePrice(BaseTool):
    __version__ = "0.1.0"
    name = "MoleculePrice"
    func_name = 'get_molecule_price'
    description = "Get the cheapest available price of a molecule."
    categories = ["Molecule"]
    tags = ["Molecular Information", "ChemSpace"]
    required_envs = [("CHEMSPACE_API_KEY", "The API key for ChemSpace.")]
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    output_sig = [('price', 'str', 'Description of the cheapest available price of the molecule')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'price': '25g of this molecule cost 143 USD and can be purchased at A2B Chem.'}},
    ]

    def __init__(self, chemspace_api_key: Optional[str] = None, init=True, interface='code') -> None:
        chemspace_api_key = os.getenv("CHEMSPACE_API_KEY", None)
        if chemspace_api_key is None:
            raise ChemMCPApiNotFoundError("CHEMSPACE_API_KEY environment variable not set.")
        self.chemspace = ChemSpace(chemspace_api_key)
        super().__init__(init=init, interface=interface)

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMCPInputError(f"smiles `{smiles}` is not a valid SMILES string.")
        price = self.chemspace.buy_mol(smiles)
        return price


if __name__ == "__main__":
    run_mcp_server()
