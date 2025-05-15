from abc import ABC, abstractmethod
import logging
import inspect
import functools
from typing import Dict, List, Tuple, Literal, Any
# from pydantic import BaseModel, Field, field_validator


logger = logging.getLogger(__name__)


# TODO: Add a validator for the required class attrs
# class _BaseToolMeta(ABC):
#     __version__: str = Field(..., description="The version of the tool.", pattern=r"^\d+\.\d+\.\d+$")
#     name: str = Field(..., description="The name of the tool.", min_length=1)
#     func_name: str = Field(..., description="The function name of the tool.", min_length=1)
#     description: str = Field(..., description="The description of the tool.", min_length=1)
#     categories: List[Literal["Molecule", "Reaction", "General"]] = Field(..., description="The categories of the tool.", min_length=1)
#     tags: List[str] = Field(..., description="The tags of the tool.", min_length=1)
#     required_envs: List[Tuple[str, str]] = Field(..., description="The required environment variables for the tool.")
#     text_input_sig: List[Tuple[str, str, str, str]] = Field(..., description="The text input signature of the tool. Each element is a tuple of (arg_name, arg_type, arg_default, arg_description).")
#     code_input_sig: List[Tuple[str, str, str, str]] = Field(..., description="The code input signature of the tool. Each element is a tuple of (arg_name, arg_type, arg_default, arg_description).")
#     output_sig: List[Tuple[str, str, str]] = Field(..., description="The output signature of the tool. Each element is a tuple of (output_name, output_type, output_description).")
#     examples: List[Dict[Literal["text_input", "code_input", "output"], Dict[str, Any]]] = Field(..., description="The examples of the tool.")


class BaseTool(ABC):
    _registered_tool = True
    _registered_mcp_tool = False
    _required_class_attrs = ['__version__', 'name', 'func_name', 'description', 'categories', 'tags', 'required_envs', 'text_input_sig', 'code_input_sig', 'output_sig', 'examples']
    
    __version__: str
    name: str
    func_name: str
    description: str
    categories: List[Literal["Molecule", "Reaction", "General"]]
    tags: List[str]
    required_envs: List[Tuple[str, str]]  # [(env_name, env_description), ...]
    code_input_sig: List[Tuple]  # [("arg_name", "arg_type", "arg_default", "arg_description"), ...]
    text_input_sig: List[Tuple]  # [("arg_name", "arg_type", "arg_default", "arg_description"), ...]
    output_sig: List[Tuple]  # [("output_name", "output_type", "output_description"), ...]
    examples: List[Dict[Literal["text_input", "code_input", "output"], Dict[str, Any]]]

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        # Check if the required class attributes are defined
        missing = [attr for attr in cls._required_class_attrs if not hasattr(cls, attr)]
        if missing:
            raise TypeError(
                f"{cls.__name__} must define class attrs: {', '.join(missing)}"
            )
        
        # Check if the version string is valid
        version = getattr(cls, "__version__", None)
        if version is None:
            raise TypeError(f"{cls.__name__} must define an valid version string as __version__ class attribute.")
        try:
            version_parts = version.split('.')
            if len(version_parts) != 3:
                raise ValueError(f"Invalid version string: {version}. Version must be in the format 'major.minor.patch'.")
            for part in version_parts:
                if not part.isdigit():
                    raise ValueError(f"Invalid version string: {version}. Version must be in the format 'major.minor.patch'.")
        except AttributeError:
            raise TypeError(f"{cls.__name__} must define an valid version string as __version__ class attribute.")

    @classmethod
    def get_doc(cls, interface='code'):
        assert interface in ('text', 'code'), "Interface '%s' is not supported. Please use 'text' or 'code'." % interface
        inputs = ""
        for name, type_, description in (cls.code_input_sig if interface == 'code' else cls.text_input_sig):
            inputs += f"    {name} ({type_}): {description}\n"
        outputs = ""
        for name, type_, description in cls.output_sig:
            outputs += f"    {name} ({type_}): {description}\n"
        
        doc = f"""{cls.description}

Args:
{inputs}
Returns:
{outputs}
"""
        return doc

    def __init__(self, init=True, interface='code') -> None:
        assert interface in ('text', 'code'), "Interface '%s' is not supported. Please use 'text' or 'code'." % interface
        self.interface = interface
        super().__init__()
        if init:
            self._init_modules()

    def _init_modules(self):
        pass

    def __call__(self, *args, **kwargs):
        logger.debug("Calling `{}` in __call__".format(self.__class__.name))
        if self.interface == 'text':
            r = self.run_text(args[0])
        elif self.interface == 'code':
            r = self.run_code(*args, **kwargs)
        else:
            raise NotImplementedError("Interface '%s' is not supported. Please use 'text' or 'code'." % self.interface)
        logger.debug("Ending `{}` in __call__".format(self.__class__.name))
        return r

    def run_text(self, query, *args, **kwargs):
        return self._run_text(query, *args, **kwargs)
    
    def _run_text(self, query, *args, **kwargs):
        # Check the signature of self._run_base
        sig = inspect.signature(self._run_base)
        params = list(sig.parameters.values())[1:]
        
        # If the number of parameters is 1, and the type is str, then we assume the tool is text-compatible
        if len(params) == 1 and params[0].annotation == str:
            return self._run_base(query)
        else:
            raise NotImplementedError("Text interface is not implemented for this tool yet.")
    
    def run_code(self, *args, **kwargs):
        return self._run_code(*args, **kwargs)
    
    def _run_code(self, *args, **kwargs):
        return self._run_base(*args, **kwargs)
    
    @abstractmethod
    def _run_base(self, *args, **kwargs):
        raise NotImplementedError
    

class ChemMCPManager:
    _tools: list[type] = []

    @staticmethod
    def register_tool(cls):
        cls._registered_mcp_tool = True
        ChemMCPManager._tools.append(cls)
        return cls
    
    @staticmethod
    def get_registered_tools():
        return list(ChemMCPManager._tools)
    
    @staticmethod
    def init_mcp_tools(mcp):
        for cls in ChemMCPManager.get_registered_tools():
            inst = cls()

            sig = inspect.signature(cls._run_base)
            params = list(sig.parameters.values())[1:]
            exposed_sig = sig.replace(parameters=params)

            def make_wrapper(inst, func):
                @functools.wraps(func)
                def wrapper(*args, **kwargs):
                    return inst.run_code(*args, **kwargs)
                return wrapper

            wrapper = make_wrapper(inst, cls._run_base)
            wrapper.__signature__ = exposed_sig
            wrapper.__name__      = cls.func_name
            wrapper.__doc__       = cls.get_doc(interface='code')

            existing = set(mcp._tool_manager._tools.keys())
            if cls.func_name not in existing:
                mcp.tool()(wrapper)
        
        logger.info(f"Initialized {len(ChemMCPManager._tools)} tools to MCP.")


# def register_mcp_tool(mcp, is_async: bool = False):
#     """
#     Class decorator to auto-register a BaseTool subclass as an MCP tool.

#     Args:
#         mcp:       The MCP instance (e.g. imported from your mcp_app).
#         is_async:  If True, creates an async wrapper (not yet implemented).
#                    If False, creates a sync wrapper.
#     """
#     if is_async:
#         raise NotImplementedError("Async registration not supported yet.")

#     def decorator(cls):
#         instance = cls()  # __init__ still runs as defined in your class

#         # 1. Inspect the signature of the public entrypoint (_run_base)
#         original_sig = inspect.signature(cls._run_base)

#         # 2. Drop the first 'self' parameter so MCP sees only actual tool args
#         params = list(original_sig.parameters.values())[1:]
#         exposed_sig = original_sig.replace(parameters=params)

#         # 3. Build a sync wrapper that instantiates & calls run_code()
#         def wrapper(*args, **kwargs):
#             return instance.run_code(*args, **kwargs)

#         # 4. Copy name, docstring, annotations, etc. from _run_base
#         wrapper = functools.wraps(cls._run_base)(wrapper)
#         wrapper.__signature__ = exposed_sig
#         wrapper.__name__      = cls.func_name
#         wrapper.__doc__       = cls.get_doc(interface='code')

#         # 5. Register with MCP if not already registered
#         existing = set(mcp._tool_manager._tools.keys())
#         if cls.func_name not in existing:
#             mcp.tool()(wrapper)

#         return cls

#     return decorator
