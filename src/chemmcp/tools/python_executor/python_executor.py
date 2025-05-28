import os
import logging
from typing import Optional

from ...utils.base_tool import BaseTool
from ...utils.mcp_app import ChemMCPManager, run_mcp_server
from .jupyter_backbone import JupyterBackbone


logger = logging.getLogger(__name__)


@ChemMCPManager.register_tool
class PythonExecutor(BaseTool):
    __version__ = "0.1.0"
    name = "PythonExecutor"
    func_name = 'run_code'
    description = "Execute Python code in a Jupyter notebook. New packages can be installed by running `!pip install <package_name>`."
    categories = ["General"]
    tags = ["Code Execution"]
    required_envs = []
    text_input_sig = [("code", "str", "N/A", "The Python code to execute.")]
    code_input_sig = [("code", "str", "N/A", "The Python code to execute.")]
    output_sig = [("result", "str", "The result of the Python code execution.")]
    examples = [
        {'text_input': {'code': 'print("Hello, world!")'}, 'code_input': {'code': 'print("Hello, world!")'}, 'output': {'result': 'Hello, world!'}},
    ]

    def __init__(self, kernel_url: Optional[str] = None, init=True, interface='code'):
        kernel_url = kernel_url or os.getenv("JUPYTER_KERNEL_GATEWAY_URL", None)
        self.jupyter = JupyterBackbone(kernel_url=kernel_url)
        super().__init__(init, interface=interface)

    def __del__(self):
        self.jupyter.close()

    def _run_base(self, code: str) -> str:
        result = self.jupyter.run_code(code)
        return result


if __name__ == "__main__":
    run_mcp_server()
