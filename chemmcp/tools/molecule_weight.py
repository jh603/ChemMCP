from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from ..utils.base_tool import BaseTool, register_mcp_tool
from ..utils.mcp_app import mcp_instance, run_mcp_server


@register_mcp_tool(mcp_instance)
class MoleculeWeight(BaseTool):
    __version__ = "0.1.0"
    name = "MoleculeWeight"
    func_name = 'cal_molecular_weight'
    description = "Calculate molecular weight."
    categories = ["Molecule"]
    tags = ["Molecular Property", "RDKit"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule')]
    output_sig = [('weight', 'float', 'Molecular weight of the molecule')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'weight': 46.041864812}},
    ]

    def _run_base(self, smiles: str) -> float:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return f"Error: `{smiles}` is not a valid SMILES string."
        mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
        return mol_weight


if __name__ == "__main__":
    run_mcp_server()

