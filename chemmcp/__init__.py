import importlib

_tool_module_map = {
    "WebSearch":    "tools.web_search",
    "MoleculeCaptioner": "molecule_description",
    "MoleculeGenerator": "molecule_description",
    "PubchemSearchQA": "pubchem_search",
    "PubchemSearch": "pubchem_search",  # Not registerred as an MCP tool
    "ForwardSynthesis": "rxn4chem",
    "Retrosynthesis": "rxn4chem",
    "MolSimilarity": "molecule_ops",
    "SMILES2Weight": "molecule_ops",
    "FuncGroups": "molecule_ops",
    "CanonicalizeSMILES": "molecule_ops",
    "CountMolAtoms": "molecule_ops",
    "GetMoleculePrice": "molecule_info",
    "PatentCheck": "molecule_info",
    "SolubilityPredictor": "property_prediction",
    "LogDPredictor": "property_prediction",
    "BBBPPredictor": "property_prediction",
    "ToxicityPredictor": "property_prediction",
    "HIVInhibitorPredictor": "property_prediction",
    "SideEffectPredictor": "property_prediction",
    "IUPAC2SMILES": "name_conversion",
    "SMILES2IUPAC": "name_conversion",
    "SMILES2Formula": "name_conversion",
    "Name2SMILES": "name_conversion",
    "SELFIES2SMILES": "name_conversion",
    "SMILES2SELFIES": "name_conversion",
}

__all__ = list(_tool_module_map.keys())

def __getattr__(name: str):
    if name in __all__:
        module_name = _tool_module_map.get(name)
        if module_name is None:
            raise AttributeError(f"No mapping for tool {name!r} in chemmtk")
        module = importlib.import_module(f"{__name__}.{module_name}")
        try:
            return getattr(module, name)
        except AttributeError:
            raise ImportError(f"Module {module_name!r} has no attribute {name!r}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return list(__all__)
