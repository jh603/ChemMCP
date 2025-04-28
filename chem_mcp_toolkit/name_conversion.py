import os
import logging
from rdkit import Chem
import rdkit.Chem.rdMolDescriptors as molD
import pubchempy as pcp
import selfies as sf

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


@mcp.tool()
async def convert_iupac_to_smiles(iupac: str) -> str:
    """Convert IUPAC name to SMILES string.

    Args:
        iupac: IUPAC name of the molecule
    Returns:
        str: SMILES string of the molecule
    """

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


@mcp.tool()
async def convert_smiles_to_iupac(smiles: str) -> str:
    """Convert SMILES to IUPAC name.

    Args:
        smiles: SMILES string of the molecule
    Returns:
        str: IUPAC name of the molecule
    """
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


@mcp.tool()
def convert_smiles_to_formula(smiles: str) -> str:
    """Convert SMILES to molecular formula.

    Args:
        smiles: SMILES string of the molecule
    Returns:
        str: Molecular formula of the molecule
    """
    if not is_smiles(smiles):
        raise ChemMTKInputError("The input is not a valid SMILES string.")
    
    try:
        formula = smiles2formula(smiles)
    except KeyboardInterrupt:
        raise
    except Exception as e:
        raise ChemMTKToolProcessError("Failed to process the SMILES string.") from e
    
    return formula


@mcp.tool()
async def convert_chemical_name_to_smiles(name: str):
    """Convert chemical name to SMILES string.

    Args:
        name: Chemical name of the molecule
    Returns:
        str: SMILES string of the molecule
    """
    try:
        smi = pubchem_name2smiles(name)
        logger.debug("Looking up PubChem succeeded.")
    except ChemMTKSearchFailError as e:
        logger.debug("Looking up PubChem failed.")
        raise e
    
    return smi


@mcp.tool()
def convert_selfies_to_smiles(selfies: str) -> str:
    """Convert SELFIES to SMILES string.

    Args:
        selfies: SELFIES string of the molecule
    Returns:
        str: SMILES string of the molecule
    """
    try:
        smiles = sf.decoder(selfies)
    except KeyboardInterrupt:
        raise
    except:
        raise ChemMTKToolProcessError("Cannot convert the SELFIES into SMILES, possibly because it is not a valid SELFIES string.")
    return smiles


@mcp.tool()
def convert_smiles_to_selfies(smiles: str) -> str:
    """Convert SMILES to SELFIES string.

    Args:
        smiles: SMILES string of the molecule
    Returns:
        str: SELFIES string of the molecule
    """
    if not is_smiles(smiles):
        raise ChemMTKInputError("The input is not a valid SMILES string.")
    try:
        selfies = sf.encoder(smiles)
    except KeyboardInterrupt:
        raise
    except:
        raise ChemMTKToolProcessError("Cannot convert the SMILES into SELFIES, possibly because it is not a valid SMILES string.")
    return selfies


# build a Starlette/uvicorn app
app = mcp.sse_app()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8001, log_level="info")
