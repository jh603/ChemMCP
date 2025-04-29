from abc import ABC, abstractmethod
import logging
import inspect
import functools
from typing import Dict, List, Tuple


logger = logging.getLogger(__name__)


class BaseTool(ABC):
    name: str
    func_name: str
    description: str
    text_input: List[Tuple]  # [("arg_name", "arg_type", "arg_description"), ...]
    code_input: List[Tuple]  # [("arg_name", "arg_type", "arg_description"), ...]
    tool_output: List[Tuple]  # [("output_name", "output_type", "output_description"), ...]
    examples: list
    
    _required_class_attrs = ['name', 'func_name', 'description', 'text_input', 'code_input', 'tool_output']

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        missing = [attr for attr in cls._required_class_attrs if not hasattr(cls, attr)]
        if missing:
            raise TypeError(
                f"{cls.__name__} must define class attrs: {', '.join(missing)}"
            )

    @classmethod
    def get_doc(cls, interface='code'):
        assert interface in ('text', 'code'), "Interface '%s' is not supported. Please use 'text' or 'code'." % interface
        inputs = ""
        for name, type_, description in (cls.code_input if interface == 'code' else cls.text_input):
            inputs += f"    {name} ({type_}): {description}\n"
        outputs = ""
        for name, type_, description in cls.tool_output:
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
        return self._run_base(query, *args, **kwargs)
    
    def run_code(self, *args, **kwargs):
        return self._run_base(*args, **kwargs)
    
    @abstractmethod
    def _run_base(self, *args, **kwargs):
        raise NotImplementedError


import inspect
import functools

def register_mcp_tool(mcp, is_async: bool = False):
    """
    Class decorator to auto-register a BaseTool subclass as an MCP tool.

    Args:
        mcp:       The MCP instance (e.g. imported from your mcp_app).
        is_async:  If True, creates an async wrapper (not yet implemented).
                   If False, creates a sync wrapper.
    """
    if is_async:
        raise NotImplementedError("Async registration not supported yet.")

    def decorator(cls):
        instance = cls()  # __init__ still runs as defined in your class

        # 1. Inspect the signature of the public entrypoint (_run_base)
        original_sig = inspect.signature(cls._run_base)

        # 2. Drop the first 'self' parameter so MCP sees only actual tool args
        params = list(original_sig.parameters.values())[1:]
        exposed_sig = original_sig.replace(parameters=params)

        # 3. Build a sync wrapper that instantiates & calls run_code()
        def wrapper(*args, **kwargs):
            return instance.run_code(*args, **kwargs)

        # 4. Copy name, docstring, annotations, etc. from _run_base
        wrapper = functools.wraps(cls._run_base)(wrapper)
        wrapper.__signature__ = exposed_sig
        wrapper.__name__      = cls.func_name
        wrapper.__doc__       = cls.get_doc(interface='code')

        # 5. Register with MCP if not already registered
        existing = set(mcp._tool_manager._tools.keys())
        if cls.func_name not in existing:
            mcp.tool()(wrapper)

        return cls

    return decorator
