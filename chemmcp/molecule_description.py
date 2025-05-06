import argparse

from transformers import T5Tokenizer, T5ForConditionalGeneration

from .utils.base_tool import BaseTool, register_mcp_tool
from .utils.smiles import is_smiles
from .utils.errors import ChemMTKInputError
from .mcp_app import mcp_instance


@register_mcp_tool(mcp_instance)
class MoleculeCaptioner(BaseTool):
    name = "MoleculeCaptioner"
    func_name = "generate_molecule_caption"
    description = "Generate a textual description of the molecule from its SMILES representation with MolT5. This tool uses neural networks to generate descriptions, which may not be accurate or correct. Please first try other tools that provide accurate and authoritative information, and only use this one as the last resort."
    text_input_sig = [('smiles', 'str', 'SMILES representation of the molecule.')]
    code_input_sig = [('smiles', 'str', 'SMILES representation of the molecule.')]
    output_sig = [('description', 'str', 'Textual description of the molecule.')]
    examples = [
        {'code_input': {'smiles': 'CCO'}, 'text_input': {'smiles': 'CCO'}, 'output': {'description': 'The molecule is an ether in which the oxygen atom is linked to two ethyl groups. It has a role as an inhalation anaesthetic, a non-polar solvent and a refrigerant. It is a volatile organic compound and an ether.\n\nNote: This is a generated description and may not be accurate. Please double check the result.'}},
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
    
    def _run_base(self, smiles: str) -> str:
        if not is_smiles(smiles):
            raise ChemMTKInputError("The input is not a valid SMILES string.")
        
        return self._run_molt5(smiles) + "\n\nNote: This is a generated description and may not be accurate. Please double check the result."


@register_mcp_tool(mcp_instance)
class MoleculeGenerator(BaseTool):
    name = "MoleculeGenerator"
    func_name = "generate_molecule_from_description"
    description = "Generate a molecule represented in SMILES with MolT5 that matches the given textual description."
    text_input_sig = [('description', 'str', 'Textual description of the molecule.')]
    code_input_sig = [('description', 'str', 'Textual description of the molecule.')]
    output_sig = [('smiles', 'str', 'SMILES representation of the molecule.')]
    examples = [
        {'code_text': {'description': 'The molecule is an ether in which the oxygen atom is linked to two ethyl groups. It has a role as an inhalation anaesthetic, a non-polar solvent and a refrigerant. It is a volatile organic compound and an ether.'}, 'text_input': {'description': 'The molecule is an ether in which the oxygen atom is linked to two ethyl groups. It has a role as an inhalation anaesthetic, a non-polar solvent and a refrigerant. It is a volatile organic compound and an ether.'}, 'output': {'smiles': 'CCO\n\nNote: This is a generated SMILES and may not be accurate. Please double check the result.'}},
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
    
    def _run_base(self, description: str) -> str:
        return self._run_molt5(description) + "\n\nNote: This is a generated SMILES and may not be accurate. Please double check the result."


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
