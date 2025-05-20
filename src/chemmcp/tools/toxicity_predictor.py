import logging

from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.property_prediction import PropertyPredictor


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class ToxicityPredictor(PropertyPredictor):
    __version__ = "0.1.0"
    name = "ToxicityPredictor"
    func_name = 'predict_toxicity'
    description = "Predict the toxicity of a molecule given its SMILES representation."
    categories = ["Molecule"]
    tags = ["Molecular Information", "Molecular Properties", "SMILES", "Neural Networks"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    output_sig = [('toxicity', 'str', 'The probability of the compound to be toxic.')]
    examples = [
        {'code_input': {'smiles': 'COC[C@@H](NC(C)=O)C(=O)NCC1=CC=CC=C1'}, 'text_input': {'smiles': 'COC[C@@H](NC(C)=O)C(=O)NCC1=CC=CC=C1'}, 'output': {'toxicity': 'The probability of the compound to be toxic is 6.04%, which means it\'s unlikely to happen.\nNote that the result is predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed.'}},
    ]

    def __init__(
        self, 
        init=True, 
        interface='code'
    ):
        super().__init__('clintox', init,  interface=interface)

    def _run_base(self, smiles: str) -> str:
        r = super()._run_base(smiles)
        return 'The probability of the compound to be toxic is {:.2f}%, which means it\'s {} to happen.'.format(r * 100, 'likely' if r >= 0.5 else 'unlikely') + '\nNote that the result is predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed.'


if __name__ == "__main__":
    run_mcp_server()
