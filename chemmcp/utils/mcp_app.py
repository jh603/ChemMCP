import argparse

from mcp.server.fastmcp import FastMCP

from .errors import catch_errors
from .base_tool import ChemMCPManager

mcp_instance = FastMCP("ChemMcpToolkit", request_timeout=300)


original_tool = mcp_instance.tool

def tool_with_catch(*dargs, **dkwargs):
    def decorator(fn):
        # first wrap errors, then register as a tool
        wrapped = catch_errors(fn)
        return original_tool(*dargs, **dkwargs)(wrapped)
    return decorator

mcp_instance.tool = tool_with_catch


def run_mcp_server():
    parser = argparse.ArgumentParser(description="Run the MCP server.")
    parser.add_argument('--sse', action='store_true', help="Run the server with SSE (Server-Sent Events) support.")
    args = parser.parse_args()

    ChemMCPManager.init_mcp_tools(mcp_instance)

    if args.sse:
        # build a Starlette/uvicorn app
        app = mcp_instance.sse_app()
        import uvicorn
        uvicorn.run(app, host="127.0.0.1", port=8001)
    else:
        # Run the MCP server with standard input/output
        mcp_instance.run(transport='stdio')
