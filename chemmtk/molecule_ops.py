import argparse

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from .utils.base_tool import BaseTool, register_mcp_tool
from .utils.errors import ChemMTKInputError
from .utils.smiles import tanimoto, is_smiles
from .utils.canonicalization import canonicalize_molecule_smiles
from .mcp_app import mcp


# Constants
# List obtained from https://github.com/rdkit/rdkit/blob/master/Data/FunctionalGroups.txt
DICT_FGS = {
    "furan": "o1cccc1",
    "aldehydes": " [CX3H1](=O)[#6]",
    "esters": " [#6][CX3](=O)[OX2H0][#6]",
    "ketones": " [#6][CX3](=O)[#6]",
    "amides": " C(=O)-N",
    "thiol groups": " [SH]",
    "alcohol groups": " [OH]",
    "methylamide": "*-[N;D2]-[C;D3](=O)-[C;D1;H3]",
    "carboxylic acids": "*-C(=O)[O;D1]",
    "carbonyl methylester": "*-C(=O)[O;D2]-[C;D1;H3]",
    "terminal aldehyde": "*-C(=O)-[C;D1]",
    "amide": "*-C(=O)-[N;D1]",
    "carbonyl methyl": "*-C(=O)-[C;D1;H3]",
    "isocyanate": "*-[N;D2]=[C;D2]=[O;D1]",
    "isothiocyanate": "*-[N;D2]=[C;D2]=[S;D1]",
    "nitro": "*-[N;D3](=[O;D1])[O;D1]",
    "nitroso": "*-[N;R0]=[O;D1]",
    "oximes": "*=[N;R0]-[O;D1]",
    "Imines": "*-[N;R0]=[C;D1;H2]",
    "terminal azo": "*-[N;D2]=[N;D2]-[C;D1;H3]",
    "hydrazines": "*-[N;D2]=[N;D1]",
    "diazo": "*-[N;D2]#[N;D1]",
    "cyano": "*-[C;D2]#[N;D1]",
    "primary sulfonamide": "*-[S;D4](=[O;D1])(=[O;D1])-[N;D1]",
    "methyl sulfonamide": "*-[N;D2]-[S;D4](=[O;D1])(=[O;D1])-[C;D1;H3]",
    "sulfonic acid": "*-[S;D4](=O)(=O)-[O;D1]",
    "methyl ester sulfonyl": "*-[S;D4](=O)(=O)-[O;D2]-[C;D1;H3]",
    "methyl sulfonyl": "*-[S;D4](=O)(=O)-[C;D1;H3]",
    "sulfonyl chloride": "*-[S;D4](=O)(=O)-[Cl]",
    "methyl sulfinyl": "*-[S;D3](=O)-[C;D1]",
    "methyl thio": "*-[S;D2]-[C;D1;H3]",
    "thiols": "*-[S;D1]",
    "thio carbonyls": "*=[S;D1]",
    "halogens": "*-[#9,#17,#35,#53]",
    "t-butyl": "*-[C;D4]([C;D1])([C;D1])-[C;D1]",
    "tri fluoromethyl": "*-[C;D4](F)(F)F",
    "acetylenes": "*-[C;D2]#[C;D1;H]",
    "cyclopropyl": "*-[C;D3]1-[C;D2]-[C;D2]1",
    "ethoxy": "*-[O;D2]-[C;D2]-[C;D1;H3]",
    "methoxy": "*-[O;D2]-[C;D1;H3]",
    "side-chain hydroxyls": "*-[O;D1]",
    "ketones": "*=[O;D1]",
    "primary amines": "*-[N;D1]",
    "nitriles": "*#[N;D1]",
}


@register_mcp_tool(mcp)
class MolSimilarity(BaseTool):
    name = "MolSimilarity"
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


@register_mcp_tool(mcp)
class SMILES2Weight(BaseTool):
    name = "SMILES2Weight"
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
    

@register_mcp_tool(mcp)
class FuncGroups(BaseTool):
    name = "FuncGroups"
    func_name = 'get_functional_groups'
    description = "Get the functional groups in a molecule."
    code_input_sig = [('smiles', 'str', 'SMILES string of the molecule.')]
    text_input_sig = [('smiles', 'str', 'SMILES string of the molecule.')]
    output_sig = [('fgs', 'str', 'A description of functional groups in the molecule.')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'fgs': 'This molecule contains alcohol groups, and side-chain hydroxyls.'}},
    ]

    def _run_base(self, smiles: str) -> str:
        def _is_fg_in_mol(mol, fg):
            fgmol = Chem.MolFromSmarts(fg)
            mol = Chem.MolFromSmiles(mol.strip())
            return len(Chem.Mol.GetSubstructMatches(mol, fgmol, uniquify=True)) > 0
        
        try:
            fgs_in_molec = [
                name
                for name, fg in DICT_FGS.items()
                if _is_fg_in_mol(smiles, fg)
            ]
            if len(fgs_in_molec) > 1:
                return f"This molecule contains {', '.join(fgs_in_molec[:-1])}, and {fgs_in_molec[-1]}."
            elif len(fgs_in_molec) == 1:
                return f"This molecule contains {fgs_in_molec[0]}."
            else:
                return "This molecule does not contain any functional groups."
        except:
            raise ChemMTKInputError("Invalid SMILES string.")


@register_mcp_tool(mcp)
class CanonicalizeSMILES(BaseTool):
    name = "CanonicalizeSMILES"
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
    

@register_mcp_tool(mcp)
class CountMolAtoms(BaseTool):
    name = "CountMolAtoms"
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
        app = mcp.sse_app()
        import uvicorn
        uvicorn.run(app, host="127.0.0.1", port=8001)
    else:
        # Run the MCP server with standard input/output
        mcp.run(transport='stdio')

