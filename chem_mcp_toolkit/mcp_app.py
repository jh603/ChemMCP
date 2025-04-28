
from mcp.server.fastmcp import FastMCP

from .utils.errors import catch_errors

mcp = FastMCP("ChemMcpToolkit", request_timeout=300)


original_tool = mcp.tool

def tool_with_catch(*dargs, **dkwargs):
    def decorator(fn):
        # first wrap errors, then register as a tool
        wrapped = catch_errors(fn)
        return original_tool(*dargs, **dkwargs)(wrapped)
    return decorator

mcp.tool = tool_with_catch
