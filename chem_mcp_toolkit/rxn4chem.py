import logging
from time import sleep
import os

from rxn4chemistry import RXN4ChemistryWrapper  # type: ignore

from .utils.base_tool import BaseTool
from .utils.errors import ChemMTKToolProcessError, ChemMTKToolInitError, ChemMTKInputError, ChemMTKApiNotFoundError
from .utils.smiles import is_smiles

from .mcp_app import mcp


logger = logging.getLogger(__name__)

class RXN4Chem(BaseTool):
    """Wrapper for RXN4Chem functionalities."""

    name: str
    description: str

    base_url: str = "https://rxn.res.ibm.com"
    sleep_time: int = 5

    rxn4chem_chemistry_wrapper = None

    def __init__(self, rxn4chem_api_key, init=True, interface='text'):
        """Init object."""
        super().__init__(init, interface=interface)

        self.rxn4chem_api_key = rxn4chem_api_key
        if RXN4Chem.rxn4chem_chemistry_wrapper is None:
            RXN4Chem.rxn4chem_chemistry_wrapper = RXN4ChemistryWrapper(
                api_key=self.rxn4chem_api_key, base_url=RXN4Chem.base_url
            )
            RXN4Chem.rxn4chem_chemistry_wrapper.create_project('ChemMTK')
        self.rxn4chem = RXN4Chem.rxn4chem_chemistry_wrapper
        if init:
            if self.rxn4chem.project_id is None:
                raise ChemMTKToolInitError("The RXN4Chem project ID cannot be initialized.")
        logger.debug("RXN4Chem project ID: %s" % self.rxn4chem.project_id)

    @staticmethod
    def retry(times: int, exceptions, sleep_time: int = 5):
        """
        Retry Decorator.

        Retries the wrapped function/method `times` times if the exceptions
        listed in ``exceptions`` are thrown
        :param times: The number of times to repeat the wrapped function/method
        :type times: Int
        :param Exceptions: Lists of exceptions that trigger a retry attempt
        :type Exceptions: Tuple of Exceptions
        """

        def decorator(func):
            def newfn(*args, **kwargs):
                attempt = 0
                while attempt < times:
                    try:
                        sleep(sleep_time)
                        return func(*args, **kwargs)
                    except exceptions:
                        print(
                            "Exception thrown when attempting to run %s, "
                            "attempt %d of %d" % (func, attempt, times)
                        )
                        attempt += 1
                return func(*args, **kwargs)

            return newfn

        return decorator


class ForwardSynthesis(RXN4Chem):
    """Predict reaction."""

    name = "ForwardSynthesis"
    func_name = "do_forward_synthesis"
    description = (
        "Predict the product of a chemical reaction. "
        "Input the SMILES of the reactants and reagents separated by a dot '.', returns SMILES of the products."
    )
    func_doc = ("reactants: str", "str")
    func_description = description
    examples = [
        {'input': 'CCN.CN1C=CC=C1C=O', 'output': 'CCNCc1cccn1C'},
    ]

    def _run_text(self, reactants: str) -> str:
        return self._run_base(reactants)

    def _run_base(self, reactants: str, *args, **kwargs) -> str:
        """Run reaction prediction."""
        # Check that input is smiles
        if not is_smiles(reactants):
            raise ChemMTKInputError("The input contains invalid SMILES. Please double-check.")
        if '.' not in reactants:
            raise ChemMTKInputError("This tool only support inputs with at least two reactants and reagents separated by a dot '.'. Please double-check.")

        try:
            prediction_id = self.predict_reaction(reactants)
            results = self.get_results(prediction_id)
            product = results["productMolecule"]["smiles"]
        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise ChemMTKToolProcessError("Failed to predict the products for the input string. Please make sure the input is a valid SMILES string containing reactants and reagents separated by a dot '.'") from e
        return product

    @RXN4Chem.retry(10, ChemMTKToolProcessError)
    def predict_reaction(self, reactants: str) -> str:
        """Make api request."""
        response = self.rxn4chem.predict_reaction(reactants)
        if "prediction_id" in response.keys():
            return response["prediction_id"]
        else:
            raise ChemMTKToolProcessError("The tool failed to predict the reaction. Maybe the input is invalid. Please make sure the input is valid SMILES of reactants separated by dot '.' and try again.")

    @RXN4Chem.retry(10, ChemMTKToolProcessError)
    def get_results(self, prediction_id: str) -> str:
        """Make api request."""
        results = self.rxn4chem.get_predict_reaction_results(prediction_id)
        if "payload" in results["response"].keys():
            return results["response"]["payload"]["attempts"][0]
        else:
            raise ChemMTKToolProcessError("Error in obtaining the results. Please make sure the input is valid SMILES of reactants separated by dot '.' and try again.")


class Retrosynthesis(RXN4Chem):
    """Predict single-step retrosynthesis."""

    name = "Retrosynthesis"
    func_name = "do_retrosynthesis"
    description = (
        "Conduct single-step retrosynthesis."
        "Input SMILES of product, returns SMILES of potential reactants separated by a dot '.' as well as the confidence. Will output multiple sets of reactants if applicable."
    )
    func_doc = ("product: str", "str")
    func_description = description
    examples = [
        {'input': 'CCO', 'output': 'There are 13 possible sets of reactants for the given product:\n1.\tReactants: C1CCOC1.CCNC(=O)c1cccn1C.[Li][AlH4]\tConfidence: 1.0\n2.\tReactants: CCN.CCO.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n3.\tReactants: CCN.CO.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n4.\tReactants: CCN.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n5.\tReactants: CCN.CCO.Cn1cccc1C=O.O.[BH4-].[Na+]\tConfidence: 1.0\n6.\tReactants: CCN.CO.Cn1cccc1C=O.O.[BH4-].[Na+]\tConfidence: 1.0\n7.\tReactants: C1CCOC1.CCN.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n8.\tReactants: CCN.Cl.Cn1cccc1C=O\tConfidence: 0.938\n9.\tReactants: CCN.Cn1cccc1C=O\tConfidence: 0.917\n10.\tReactants: CCN.Cl.Cn1cccc1C=O\tConfidence: 0.841\n11.\tReactants: C1CCOC1.CCN.Cn1cccc1C=O\tConfidence: 0.797\n12.\tReactants: C1CCOC1.CCN.CO.Cn1cccc1C=O\tConfidence: 0.647\n13.\tReactants: C1CCOC1.CC(=O)NCc1cccn1C.[Li][AlH4]\tConfidence: 1.0\n'},  
    ]

    def _run_base(self, target: str, *args, **kwargs) -> str:
        """Run retrosynthesis prediction."""
        # Check that input is smiles
        if not is_smiles(target):
            raise ChemMTKInputError("The input contains invalid SMILES. Please double-check.")

        try:
            prediction_id = self.predict_retrosynthesis(target)
            paths = self.get_paths(prediction_id)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            raise ChemMTKToolProcessError("Failed to predict the reactants for the input string. Please make sure the input is a valid SMILES string containing one chemical.") from e

        result = "There %s %d possible sets of reactants for the given product:\n" % (
            "are" if len(paths) > 1 else "is",
            len(paths),
        )
        result_list = []
        for idx, path in enumerate(paths, start=1):
            children_smiles, confidence = self._get_children_smiles_and_confidence(path)
            result_list.append((children_smiles, confidence))
        result_list.sort(key=lambda x: x[1], reverse=True)
        for idx, (children_smiles, confidence) in enumerate(result_list, start=1):
            result += f"{idx}.\tReactants: {children_smiles}\tConfidence: {confidence}\n"
        return result

    @RXN4Chem.retry(10, KeyError)
    def predict_retrosynthesis(self, target: str) -> str:
        """Make api request."""
        response = self.rxn4chem.predict_automatic_retrosynthesis(
            product=target,
            max_steps=1,
        )
        if "prediction_id" in response.keys():
            return response["prediction_id"]
        raise KeyError

    @RXN4Chem.retry(20, ChemMTKToolProcessError)
    def get_paths(self, prediction_id: str) -> str:
        """Make api request."""
        results = self.rxn4chem.get_predict_automatic_retrosynthesis_results(
            prediction_id
        )
        if "retrosynthetic_paths" not in results.keys():
            raise ChemMTKToolProcessError("Error in obtaining the results. Please make sure the input is valid SMILES and try again.")
        paths = results["retrosynthetic_paths"]
        if paths is not None:
            if len(paths) > 0:
                return paths
        if results["status"] == "PROCESSING":
            sleep(self.sleep_time * 2)
        raise ChemMTKToolProcessError("Error in obtaining the results. Please make sure the input is valid SMILES and try again.")
    
    def _get_children_smiles_and_confidence(self, path):
        children = path['children']
        children_smiles = []
        for child in children:
            smiles = child['smiles']
            children_smiles.append(smiles)
        return '.'.join(children_smiles), path['confidence']


forward_synthesis = None
retrosynthesis = None


@mcp.tool()
def do_forward_synthesis(reactants_and_reagents_smiles: str) -> str:
    """Predict the product of a chemical reaction. Input the SMILES of the reactants and reagents separated by a dot '.', returns SMILES of the products.
    
    Args:
        reactants_and_reagents_smiles: The SMILES of the reactants and reagents separated by a dot '.'.
    Returns:
        str: The SMILES of the products.
    """
    rxn4chem_api_key = os.environ.get("RXN4CHEM_API_KEY")
    if rxn4chem_api_key is None:
        raise ChemMTKApiNotFoundError("The RXN4Chem API key is not set. Please set the RXN4CHEM_API_KEY environment variable.")

    global forward_synthesis
    if forward_synthesis is None:
        forward_synthesis = ForwardSynthesis(rxn4chem_api_key=rxn4chem_api_key, init=True)
    return forward_synthesis(reactants_and_reagents_smiles)


@mcp.tool()
def do_retrosynthesis(product_smiles: str) -> str:
    """Conduct single-step retrosynthesis. Input SMILES of product, returns SMILES of potential reactants separated by a dot '.' as well as the confidence. Will output multiple sets of reactants if applicable.
    
    Args:
        product_smiles: The SMILES of the product.
    Returns:
        str: The SMILES of the reactants and the confidence.
    """
    rxn4chem_api_key = os.environ.get("RXN4CHEM_API_KEY")
    if rxn4chem_api_key is None:
        raise ChemMTKApiNotFoundError("The RXN4Chem API key is not set. Please set the RXN4CHEM_API_KEY environment variable.")

    global retrosynthesis
    if retrosynthesis is None:
        retrosynthesis = Retrosynthesis(rxn4chem_api_key=rxn4chem_api_key, init=True)
    return retrosynthesis(product_smiles)


# build a Starlette/uvicorn app
app = mcp.sse_app()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8001, log_level="info")
