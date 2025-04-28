from rdkit import Chem
import os

import molbloom

from .utils.chemspace import ChemSpace
from .mcp_app import mcp


@mcp.tool()
async def get_molecule_price(smiles: str) -> str:
    """Get the cheapest available price of a molecule.

    Args:
        smiles: SMILES string of the molecule
    Returns:
        str: Discription of the cheapest available price of the molecule
    """

    chemspace_api_key = os.getenv("CHEMSPACE_API_KEY", None)
    if chemspace_api_key is None:
        return "Error: CHEMSPACE_API_KEY environment variable not set."
    chemspace = ChemSpace(chemspace_api_key)
    price = chemspace.buy_mol(smiles)

    return price


@mcp.tool()
async def check_molecule_if_patented(smiles: str) -> str:
    """Get whether a molecule is patented or not.

    Args:
        smiles: SMILES string of the molecule
    Returns:
        str: "patented" or "not patented"
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Error: Invalid SMILES string."

    try:
        r = molbloom.buy(smiles, canonicalize=True, catalog="surechembl")
        if r:
            output = "patented"
        else:
            output = "not patented"
    except KeyboardInterrupt:
        raise
    return output


# build a Starlette/uvicorn app
app = mcp.sse_app()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8001, log_level="info")
