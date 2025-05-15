import importlib

_tool_module_map = {
    "WebSearch":    "tools.web_search",
    "MoleculeCaptioner": "tools.molecule_captioner",
    "MoleculeGenerator": "tools.molecule_generator",
    "PubchemSearchQA": "tools.pubchem_search_qa",
    "PubchemSearch": "tools.pubchem_search",  # Not registerred as an MCP tool
    "ForwardSynthesis": "tools.forward_synthesis",
    "Retrosynthesis": "tools.retrosynthesis",
    "MoleculeSimilarity": "tools.molecule_similarity",
    "MoleculeWeight": "tools.molecule_weight",
    "FunctionalGroups": "tools.functional_groups",
    "SmilesCanonicalization": "tools.smiles_canonicalization",
    "MoleculeAtomCount": "tools.molecule_atom_count",
    "MoleculePrice": "tools.molecule_price",
    "PatentCheck": "tools.patent_check",
    "SolubilityPredictor": "tools.solubility_predictor",
    "LogDPredictor": "tools.logd_predictor",
    "BbbpPredictor": "tools.bbbp_predictor",
    "ToxicityPredictor": "tools.toxicity_predictor",
    "HivInhibitorPredictor": "tools.hiv_inhibitor_predictor",
    "SideEffectPredictor": "tools.side_effect_predictor",
    "Iupac2Smiles": "tools.iupac2smiles",
    "Smiles2Iupac": "tools.smiles2iupac",
    "Smiles2Formula": "tools.smiles2formula",
    "Name2Smiles": "tools.name2smiles",
    "Selfies2Smiles": "tools.selfies2smiles",
    "Smiles2Selfies": "tools.smiles2selfies",
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
