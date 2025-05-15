import logging

from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.property_prediction import PropertyPredictor


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class SolubilityPredictor(PropertyPredictor):
    __version__ = "0.1.0"
    name = "SolubilityPredictor"
    func_name = 'predict_solubility'
    description = "Predict the solubility of a molecule given its SMILES representation."
    categories = ["Molecule"]
    tags = ["Molecular Properties", "Neural Networks"]
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    output_sig = [('solubility', 'float', 'Log solubility in mol/L.')]
    examples = [
        {'code_input': {'smiles': 'CC(C)Cl'}, 'text_input': {'smiles': 'CC(C)Cl'}, 'output': {'solubility': 'The log solubility in mol/L is -1.410.\nNote that the result is predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed.'}},
    ]

    def __init__(
        self, 
        init=True, 
        interface='code'
    ):
        super().__init__('esol', init, interface=interface)

    def _run_base(self, smiles: str) -> str:
        r = super()._run_base(smiles)
        return 'The log solubility in mol/L is {:.3f}.'.format(r) + '\nNote that the result is predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed.'


if __name__ == "__main__":
    run_mcp_server()
