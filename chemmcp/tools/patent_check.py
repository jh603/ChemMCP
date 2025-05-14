import os
import argparse

import molbloom

from ..utils.base_tool import BaseTool, register_mcp_tool
from ..utils.errors import ChemMTKApiNotFoundError, ChemMTKInputError
from ..utils.chemspace import ChemSpace
from ..utils.smiles import is_smiles
from ..utils.mcp_app import mcp_instance
    

@register_mcp_tool(mcp_instance)
class PatentCheck(BaseTool):
    __version__ = "0.1.0"
    name = "PatentCheck"
    func_name = 'check_molecule_if_patented'
    description = "Get whether a molecule is patented or not."
    code_input_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    output_sig = [('patent_status', 'str', '"patented" or "not patented"')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'patent_status': 'not patented'}},
    ]

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMTKInputError(f"smiles `{smiles}` is not a valid SMILES string.")

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
