import os
import logging
import argparse

from rdkit import Chem
import rdkit.Chem.rdMolDescriptors as molD
import pubchempy as pcp
import selfies as sf

from .utils.base_tool import BaseTool, register_mcp_tool
from .utils.errors import ChemMTKInputError, ChemMTKSearchFailError, ChemMTKToolProcessError, ChemMTKApiNotFoundError
from .utils.smiles import is_smiles
from .utils.pubchem import pubchem_iupac2cid, pubchem_name2cid
from .utils.chemspace import ChemSpace
from .mcp_app import mcp


logger = logging.getLogger(__name__)


def pubchem_iupac2smiles(
    query: str,
    strict: bool = False,
) -> str:
    cid = pubchem_iupac2cid(query, strict=strict)  # May be a tuple in the case of multiple components in the IUPAC name
    if not isinstance(cid, tuple):
        cid = (cid,)
    smiles = []
    for single_cid in cid:
        c = pcp.Compound.from_cid(single_cid)
        r = c.isomeric_smiles
        smiles.append(r.strip())
    r = '.'.join(smiles)

    return r


def pubchem_name2smiles(
    query: str,
) -> str:
    cid = pubchem_name2cid(query)
    c = pcp.Compound.from_cid(cid)
    r = c.isomeric_smiles

    return r


def pubchem_smiles2iupac(smi):
    """This function queries the given molecule smiles and returns iupac"""

    c = pcp.get_compounds(smi, 'smiles')
    
    if len(c) == 0:
        parts = smi.split('.')
        if len(parts) > 1:
            parts_iupac = []
            parts_cannot_find = []
            for part in parts:
                try:
                    iupac = pubchem_smiles2iupac(part)
                except ChemMTKSearchFailError:
                    parts_cannot_find.append(part)
                else:
                    parts_iupac.append(iupac)
            if len(parts_cannot_find) > 0:
                raise ChemMTKSearchFailError("Cannot find a matched molecule/compound for the following parts of the input SMILES: %s" % ', '.join(parts_cannot_find))
            else:
                r = ';'.join(parts_iupac)
        else:
            raise ChemMTKSearchFailError("Cannot find a matched molecule/compound. Please check the input SMILES.")
    elif len(c) >= 1:
        if len(c) > 1:
            logger.info("There are more than one molecules/compounds that match the input SMILES. Using the first matched one.")
        c = c[0]
    
        r = c.iupac_name
    
    if r is None or r == 'None':
        raise ChemMTKSearchFailError(f"The PubChem entry (CID: {c.cid}) does not have a valid IUPAC name recorded.")

    return r


def addHs(mol):
    mol = Chem.rdmolops.AddHs(mol, explicitOnly=True)
    return mol


def smiles2formula(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = addHs(mol)
    formula = molD.CalcMolFormula(mol)
    return formula


@register_mcp_tool(mcp)
class IUPAC2SMILES(BaseTool):
    name = "IUPAC2SMILES"
    func_name = 'convert_iupac_to_smiles'
    description = "Convert IUPAC name to SMILES string."
    code_input_sig = [('iupac', 'str', 'IUPAC name of the molecule')]
    text_input_sig = [('iupac', 'str', 'IUPAC name of the molecule')]
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


@register_mcp_tool(mcp)
class SMILES2IUPAC(BaseTool):
    name = "SMILES2IUPAC"
    func_name = 'convert_smiles_to_iupac'
    description = "Convert SMILES to IUPAC name."
    code_input_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    output_sig = [('iupac', 'str', 'IUPAC name of the molecule')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'iupac': 'ethanol'}},
    ]

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMTKInputError("The input is not a valid SMILES string.")
        
        try:
            name = pubchem_smiles2iupac(smiles)
            logger.debug("Looking up PubChem succeeded.")
        except KeyboardInterrupt:
            raise
        except ChemMTKSearchFailError as e:
            logger.debug("Looking up PubChem failed.")
            raise e
        
        # If PubChem fails, try STOUT
        # TODO: The STOUT package is not available anymore. Should be fixed later.
        
        return name


@register_mcp_tool(mcp)
class SMILES2Formula(BaseTool):
    name = "SMILES2Formula"
    func_name = 'convert_smiles_to_formula'
    description = "Convert SMILES to molecular formula."
    code_input_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    output_sig = [('formula', 'str', 'Molecular formula of the molecule')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'formula': 'C2H6O'}},
    ]

    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMTKInputError("The input is not a valid SMILES string.")
        
        try:
            formula = smiles2formula(smiles)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise ChemMTKToolProcessError("Failed to process the SMILES string.") from e
        
        return formula


@register_mcp_tool(mcp)
class Name2SMILES(BaseTool):
    name = "Name2SMILES"
    func_name = 'convert_chemical_name_to_smiles'
    description = "Convert chemical name to SMILES string."
    code_input_sig = [('name', 'str', 'Chemical name of the molecule')]
    text_input_sig = [('name', 'str', 'Chemical name of the molecule')]
    output_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    examples = [
        {'code_input': {'name': 'aspirin'}, 'text_input': {'name': 'aspirin'}, 'output': {'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'}},
    ]

    def _run_base(self, name: str) -> str:
        try:
            smi = pubchem_name2smiles(name)
            logger.debug("Looking up PubChem succeeded.")
        except ChemMTKSearchFailError as e:
            logger.debug("Looking up PubChem failed.")
            raise e
        
        return smi


@register_mcp_tool(mcp)
class SELFIES2SMILES(BaseTool):
    name = "SELFIES2SMILES"
    func_name = 'convert_selfies_to_smiles'
    description = "Convert SELFIES to SMILES string."
    code_input_sig = [('selfies', 'str', 'SELFIES string of the molecule')]
    text_input_sig = [('selfies', 'str', 'SELFIES string of the molecule')]
    output_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    examples = [
        {'code_input': {'selfies': '[C][C][O]'}, 'text_input': {'selfies': '[C][C][O]'}, 'output': {'smiles': 'CCO'}},
    ]

    def _run_base(self, selfies: str) -> str:
        try:
            smiles = sf.decoder(selfies)
        except KeyboardInterrupt:
            raise
        except:
            raise ChemMTKToolProcessError("Cannot convert the SELFIES into SMILES, possibly because it is not a valid SELFIES string.")
        return smiles


@register_mcp_tool(mcp)
class SMILES2SELFIES(BaseTool):
    name = "SMILES2SELFIES"
    func_name = 'convert_smiles_to_selfies'
    description = "Convert SMILES to SELFIES string."
    code_input_sig = [('smiles', 'str', 'SMILES string of the molecule')]
    text_input_sig = [('smiles', 'str', 'SMILES string of the molecule')]
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
