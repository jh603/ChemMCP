# chemmtk/__init__.py
import importlib

_tool_module_map = {
    "WebSearch":    "general_search",
    "MoleculeCaptioner": "molecule_description",
    "MoleculeGenerator": "molecule_description",
    "PubchemSearchQA": "pubchem_search",
    "PubchemSearch": "pubchem_search",  # Not registerred as an MCP tool
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
