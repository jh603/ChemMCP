import argparse

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from ..utils.base_tool import BaseTool, register_mcp_tool
from ..utils.mcp_app import mcp_instance


@register_mcp_tool(mcp_instance)
class MoleculeWeight(BaseTool):
    __version__ = "0.1.0"
    name = "MoleculeWeight"
    func_name = 'cal_molecular_weight'
    description = "Calculate molecular weight."
    code_input_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'SMILES string of the molecule')]
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

