import os
import argparse

import molbloom

from ..utils.base_tool import BaseTool
from ..utils.errors import ChemMCPApiNotFoundError, ChemMCPInputError
from ..tool_utils.smiles import is_smiles
from ..utils.mcp_app import ChemMCPManager, run_mcp_server
    

@ChemMCPManager.register_tool
class PatentCheck(BaseTool):
    __version__ = "0.1.0"
    name = "PatentCheck"
    func_name = 'check_molecule_if_patented'
    description = "Get whether a molecule is patented or not."
    categories = ["Molecule"]
    tags = ["Patent", "MolBloom", "Web"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    output_sig = [('patent_status', 'str', '"patented" or "not patented"')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'patent_status': 'not patented'}},
    ]

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMCPInputError(f"smiles `{smiles}` is not a valid SMILES string.")

        try:
            r = molbloom.buy(smiles, canonicalize=True, catalog="surechembl")
            if r:
                output = "patented"
            else:
                output = "not patented"
        except KeyboardInterrupt:
            raise
        return output


if __name__ == "__main__":
    run_mcp_server()
