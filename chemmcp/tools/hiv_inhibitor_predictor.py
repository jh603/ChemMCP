import logging

from ..utils.mcp_app import ChemMCPManager, run_mcp_server
from ..tool_utils.property_prediction import PropertyPredictor


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class HivInhibitorPredictor(PropertyPredictor):
    __version__ = "0.1.0"
    name = "HivInhibitorPredictor"
    func_name = 'predict_hiv'
    description = "Predict the HIV inhibition of a molecule given its SMILES representation."
    categories = ["Molecule"]
    tags = ["Molecular Properties", "Neural Networks"]
    required_envs = []
    code_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    text_input_sig = [('smiles', 'str', 'N/A', 'SMILES string of the molecule.')]
    output_sig = [('hiv', 'str', 'The probability of the compound to be an inhibitor of HIV replication.')]
    examples = [
        {'code_input': {'smiles': 'CC1=CN(C2C=CCCC2O)C(=O)NC1=O'}, 'text_input': {'smiles': 'CC1=CN(C2C=CCCC2O)C(=O)NC1=O'}, 'output': {'hiv': 'The probability of the compound to be an inhibitor of HIV replication is 6.01%, which means it\'s unlikely to happen.\nNote that the result is predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed.'}},
    ]

    def __init__(
        self, 
        init=True, 
        interface='text'
    ):
        super().__init__('hiv', init, interface=interface)

    def _run_base(self, smiles: str) -> str:
        r = super()._run_base(smiles)
        return 'The probability of the compound to be an inhibitor of HIV replication is {:.2f}%, which means it\'s {} to happen.'.format(r * 100, 'likely' if r >= 0.5 else 'unlikely') + '\nNote that the result is predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed.'


if __name__ == "__main__":
    run_mcp_server()
