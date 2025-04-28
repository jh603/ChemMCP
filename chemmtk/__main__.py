import argparse

from . import general_search
from . import molecule_ops
from . import molecule_info
from . import name_conversion
from . import molecule_description
from . import property_prediction
from . import pubchem_search
from . import rxn4chem

from .mcp_app import mcp


# if __name__ == "__main__":
#     # Initialize and run the server
#     mcp.run(transport='stdio')




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the MCP server.")
    parser.add_argument('--sse', action='store_true', help="Run the server with SSE (Server-Sent Events) support.")
    args = parser.parse_args()

    if args.sse:
        # build a Starlette/uvicorn app
        app = mcp.sse_app()
        import uvicorn
        uvicorn.run(app, host="127.0.0.1", port=8001, log_level="info")
    else:
        # Run the MCP server with standard input/output
        mcp.run(transport='stdio')
