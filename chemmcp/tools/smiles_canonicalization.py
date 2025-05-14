import argparse

from ..utils.base_tool import BaseTool, register_mcp_tool
from ..utils.errors import ChemMTKInputError
from ..utils.canonicalization import canonicalize_molecule_smiles
from ..utils.mcp_app import mcp_instance


@register_mcp_tool(mcp_instance)
class SMILESCanonicalization(BaseTool):
    __version__ = "0.1.0"
    name = "SMILESCanonicalization"
    func_name = 'canonicalize_smiles'
    description = "Canonicalize a molecular SMILES string."
    code_input_sig = [('smiles', 'str', 'SMILES string of the molecule.'), ('isomeric', 'bool', 'Whether to include isomeric information. Default is True.'), ('kekulization', 'bool', 'Whether to perform kekulization. Default is True.'), ('keep_atom_map', 'bool', 'Whether to keep atom mapping numbers, if any. Default is True.')]
    text_input_sig = [('smiles', 'str', 'SMILES string of the molecule.')]  # TODO: Support options.
    output_sig = [('canonical_smiles', 'str', 'Canonicalized SMILES string.')]
    examples = [
        {'code_input': {'smiles': 'C(O)C', 'isomeric': True, 'kekulization': True, 'keep_atom_map': False}, 'text_input': {'smiles': 'C(O)C'}, 'output': {'canonical_smiles': 'CCO'}},
    ]

    def _run_base(self, smiles: str, isomeric: bool = True, kekulization: bool = True, keep_atom_map: bool = True) -> str:
        smiles = canonicalize_molecule_smiles(smiles, isomeric=isomeric, kekulization=kekulization, keep_atom_map=keep_atom_map)
        if smiles is None:
            raise ChemMTKInputError("Invalid SMILES string.")
        return smiles


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

