import asyncio
import logging

from transformers import T5Tokenizer, T5ForConditionalGeneration

from .utils.base_tool import BaseTool
from .utils.smiles import is_smiles
from .utils.errors import ChemMTKInputError

from .mcp_app import mcp


class MoleculeCaptioner(BaseTool):
    name = "MoleculeCaptioner"
    func_name = "generate_molecule_caption"
    description = "Input the SMILES of a molecule/compound, returns the textual description of the molecule/compound. This tool uses neural networks to generate descriptions, which may not be accurate or correct. You should try the PubchemSearch or WebSearch tool first, which provide accurate and authoritative information, and only use this one when other tools cannot provides useful information."
    func_doc = ("smiles: str", "str")
    func_description = description
    examples = [
        {'input': 'CCO', 'output': 'The molecule is an ether in which the oxygen atom is linked to two ethyl groups. It has a role as an inhalation anaesthetic, a non-polar solvent and a refrigerant. It is a volatile organic compound and an ether.'},
    ]

    def __init__(self, init=True, interface='text') -> None:
        self.tokenizer, self.model = None, None
        super().__init__(init, interface=interface)

    def _init_modules(self):
        self.tokenizer, self.model = self.__load_molt5()
        
    def __load_molt5(self):
        tokenizer = T5Tokenizer.from_pretrained("laituan245/molt5-large-smiles2caption", model_max_length=1024)
        model = T5ForConditionalGeneration.from_pretrained('laituan245/molt5-large-smiles2caption')
        return tokenizer, model
    
    def _run_molt5(self, smiles):
        if self.tokenizer is None or self.model is None:
            self._init_modules()
        input_ids = self.tokenizer(smiles, return_tensors="pt").input_ids
        outputs = self.model.generate(input_ids, num_beams=5, max_length=1024)
        text = self.tokenizer.decode(outputs[0], skip_special_tokens=True)
        return text
    
    def _run_base(self, smiles, *args, **kwargs):
        if not is_smiles(smiles):
            raise ChemMTKInputError("The input is not a valid SMILES string.")
        
        return self._run_molt5(smiles) + "\n\nNote: This is a generated description and may not be accurate. Please double check the result."


class MoleculeGenerator(BaseTool):
    name = "MoleculeGenerator"
    func_name = "generate_molecule_from_description"
    description = "Input a description of a molecule/compound, returns the SMILES representation of the molecule/compound. This tool uses neural networks to generate molecules, which may not be accurate or correct."
    func_doc = ("description: str", "str")
    func_description = description
    examples = [
        {'input': 'The molecule is an ether in which the oxygen atom is linked to two ethyl groups. It has a role as an inhalation anaesthetic, a non-polar solvent and a refrigerant. It is a volatile organic compound and an ether.', 'output': 'CCO'},
    ]

    def __init__(self, init=True, interface='text') -> None:
        self.tokenizer, self.model = None, None
        super().__init__(init, interface=interface)

    def _init_modules(self):
        self.tokenizer, self.model = self.__load_molt5()
        
    def __load_molt5(self):
        tokenizer = T5Tokenizer.from_pretrained("laituan245/molt5-large-caption2smiles", model_max_length=512)
        model = T5ForConditionalGeneration.from_pretrained('laituan245/molt5-large-caption2smiles')
        return tokenizer, model
    
    def _run_molt5(self, text):
        if self.tokenizer is None or self.model is None:
            self._init_modules()
        input_ids = self.tokenizer(text, return_tensors="pt").input_ids
        outputs = self.model.generate(input_ids, num_beams=5, max_length=512)
        smiles = self.tokenizer.decode(outputs[0], skip_special_tokens=True)
        return smiles
    
    def _run_base(self, description, *args, **kwargs):
        return self._run_molt5(description) + "\n\nNote: This is a generated SMILES and may not be accurate. Please double check the result."


molecule_captioner = MoleculeCaptioner()
molecule_generator = MoleculeGenerator()

@mcp.tool()
async def generate_molecule_caption(smiles: str) -> str:
    """Generate a textual description of the molecule from its SMILES representation with MolT5.

    Args:
        smiles: SMILES representation of the molecule
    Returns:
        str: Textual description of the molecule
    """
    global molecule_captioner
    if molecule_captioner is None:
        molecule_captioner = MoleculeCaptioner()
    return molecule_captioner(smiles)

# @mcp.tool()
# async def generate_molecule_caption(smiles: str) -> str:
#     logging.debug(f"🖋️  generate_molecule_caption start: {smiles}")
#     try:
#         # your actual call
#         text = molecule_captioner(smiles)
#         logging.debug("🖋️  caption done")
#         return text
#     except Exception as e:
#         logging.exception("💥 captioner raised")
#         return f"ERROR: {e}"


@mcp.tool()
async def generate_molecule_from_description(description: str) -> str:
    """Generate a SMILES representation of the molecule from its textual description with MolT5.

    Args:
        description: Textual description of the molecule
    Returns:
        str: SMILES representation of the molecule
    """
    global molecule_generator
    if molecule_generator is None:
        molecule_generator = MoleculeGenerator()
    return molecule_generator(description)


# build a Starlette/uvicorn app
app = mcp.sse_app()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8001, log_level="info")

