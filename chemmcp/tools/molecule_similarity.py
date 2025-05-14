import argparse

from ..utils.base_tool import BaseTool, register_mcp_tool
from ..utils.errors import ChemMTKInputError
from ..utils.smiles import tanimoto, is_smiles
from ..utils.mcp_app import mcp_instance


@register_mcp_tool(mcp_instance)
class MoleculeSimilarity(BaseTool):
    __version__ = "0.1.0"
    name = "MoleculeSimilarity"
    func_name = 'cal_molecule_similarity'
    description = "Get the Tanimoto similarity of two molecules. It Can also be used to check if two molecules are identical."
    code_input_sig = [('smiles1', 'str', 'SMILES string of the first molecule'), ('smiles2', 'str', 'SMILES string of the second molecule.')]
    text_input_sig = [('smiles_pair', 'str', 'SMILES strings of the two molecules, separated by a semicolon.')]
    output_sig = [('similarity', 'str', 'Tanimoto similarity score and similarity description')]
    examples = [
        {'code_input': {'smiles1': 'CCO', 'smiles2': 'CCN'}, 'text_input': {'smiles_pair': 'CCO;CCN'}, 'output': {'similarity': 'The Tanimoto similarity between CCO and CCN is 0.3333, indicating that the two molecules are not similar.'}},
        {'code_input': {'smiles1': 'CCO', 'smiles2': 'C(O)C'}, 'text_input': 'CCO;C(O)C', 'output': {'similarity': 'Input Molecules Are Identical'}},
    ]

    def _run_base(self, smiles1: str, smiles2: str) -> str:
        if not is_smiles(smiles1):
            raise ChemMTKInputError(f"smiles1 `{smiles1}` is not a valid SMILES string.")
        
        if not is_smiles(smiles2):
            raise ChemMTKInputError(f"smiles2 `{smiles2}` is not a valid SMILES string.")

        similarity = tanimoto(smiles1, smiles2)

        if isinstance(similarity, str):
            return similarity

        sim_score = {
            0.9: "very similar",
            0.8: "similar",
            0.7: "somewhat similar",
            0.6: "not very similar",
            0: "not similar",
        }
        if similarity == 1:
            return "The input molecules are identical."
        else:
            val = sim_score[
                max(key for key in sim_score.keys() if key <= round(similarity, 1))
            ]
            message = f"The Tanimoto similarity between {smiles1} and {smiles2} is {round(similarity, 4)}, indicating that the two molecules are {val}."
        return message


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

