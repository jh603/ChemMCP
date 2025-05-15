import os
import logging

from ..utils.base_tool import BaseTool, register_mcp_tool
from ..utils.errors import ChemMTKSearchFailError, ChemMTKApiNotFoundError
from ..utils.names import pubchem_iupac2smiles
from ..tool_utils.chemspace import ChemSpace
from ..utils.mcp_app import mcp_instance, run_mcp_server


logger = logging.getLogger(__name__)


@register_mcp_tool(mcp_instance)
class Iupac2Smiles(BaseTool):
    __version__ = "0.1.0"
    name = "Iupac2Smiles"
    func_name = 'convert_iupac_to_smiles'
    description = "Convert IUPAC name to SMILES string."
    categories = ["Molecule"]
    tags = ["SMILES", "RDKit", "PubChem", "ChemSpace", "IUPAC"]
    required_envs = []
    code_input_sig = [('iupac', 'str', 'N/A', 'IUPAC name of the molecule')]
    text_input_sig = [('iupac', 'str', 'N/A', 'IUPAC name of the molecule')]
    output_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    examples = [
        {'code_input': {'iupac': 'ethanol'}, 'text_input': {'iupac': 'ethanol'}, 'output': {'smiles': 'CCO'}},
    ]

    def _run_base(self, iupac: str) -> str:
        smi = None

        try:
            smi = pubchem_iupac2smiles(iupac, strict=True)
            logger.debug("Looking up PubChem succeeded.")
            return smi
        except ChemMTKSearchFailError as e:
            logger.debug("Looking up PubChem failed.")

        # If PubChem fails, try ChemSpace
        chemspace_api_key = os.getenv("CHEMSPACE_API_KEY", None)
        if not chemspace_api_key:
            logger.debug("Looking up ChemSpace failed, because ChemSpace API is not set.")
            raise ChemMTKApiNotFoundError("Cannot find the API key for ChemSpace. Please set the CHEMSPACE_API_KEY environment variable.")
        
        chemspace = ChemSpace(chemspace_api_key)
        tmp = chemspace.convert_mol_rep(iupac, "smiles")
        try:
            smi = tmp.split(":")[1].strip()
            logger.debug("Looking up ChemSpace succeeded.")
            return smi
        except IndexError as e:
            logger.debug("Looking up ChemSpace failed, due to IndexError.")
            smi = None

        if smi is None:
            raise ChemMTKSearchFailError('Cannot find a matched molecule/compound for the input IUPAC name from PubChem or ChemSpace. This may be because the input IUPAC name is not valid or the molecule is not in the databases. Please double check the input or try other tool.') from e
        
        # If ChemSpace fails, try STOUT
        # TODO: The STOUT package is not available anymore. Should be fixed later.

        return smi


if __name__ == "__main__":
    run_mcp_server()
