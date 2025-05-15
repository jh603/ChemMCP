import selfies as sf

from ..utils.base_tool import BaseTool, register_mcp_tool
from ..utils.errors import ChemMTKInputError, ChemMTKToolProcessError
from ..utils.smiles import is_smiles
from ..utils.mcp_app import mcp_instance, run_mcp_server


@register_mcp_tool(mcp_instance)
class Smiles2Selfies(BaseTool):
    __version__ = "0.1.0"
    name = "Smiles2Selfies"
    func_name = 'convert_smiles_to_selfies'
    description = "Convert SMILES to SELFIES string."
    categories = ["Molecule"]
    tags = ["SMILES", "SELFIES"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    output_sig = [('selfies', 'str', 'SELFIES string of the molecule')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'selfies': '[C][C][O]'}},
    ]

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMTKInputError("The input is not a valid SMILES string.")
        try:
            selfies = sf.encoder(smiles)
        except KeyboardInterrupt:
            raise
        except:
            raise ChemMTKToolProcessError("Cannot convert the SMILES into SELFIES, possibly because it is not a valid SMILES string.")
        return selfies


if __name__ == "__main__":
    run_mcp_server()
