import argparse

from rdkit import Chem

from ..utils.base_tool import BaseTool, register_mcp_tool
from ..utils.errors import ChemMTKInputError
from ..utils.mcp_app import mcp_instance


@register_mcp_tool(mcp_instance)
class MoleculeAtomCount(BaseTool):
    __version__ = "0.1.0"
    name = "MolAtomCount"
    func_name = 'count_molecule_atoms'
    description = "Count the number of atoms of each type in a molecule."
    code_input_sig = [('smiles', 'str', 'SMILES string of the molecule.')]
    text_input_sig = [('smiles', 'str', 'SMILES string of the molecule.')]
    output_sig = [('atom_counts', 'str', 'A description of atom numbers in the molecule.')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'atom_counts': 'There are altogether 3 atoms (omitting hydrogen atoms). The types and corresponding numbers are: {\'C\': 2, \'O\': 1}'}},
    ]

    def _run_base(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ChemMTKInputError(f"`{smiles}` is not a valid SMILES string.")
        
        num_atoms = mol.GetNumAtoms()
        # Get the atom types
        atom_types = [atom.GetSymbol() for atom in mol.GetAtoms()]
        # Count the occurrences of each atom type
        atom_type_counts = {atom: atom_types.count(atom) for atom in set(atom_types)}
        
        text = "There are altogether %d atoms (omitting hydrogen atoms). The types and corresponding numbers are: %s" % (num_atoms, str(atom_type_counts))

        return text


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

