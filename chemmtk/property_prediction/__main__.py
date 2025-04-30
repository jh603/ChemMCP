import argparse
import logging

logger = logging.getLogger(__name__)


# Import all the modules to register them with the MCP server
from . import property_prediction

from ..mcp_app import mcp

logger.info("ChemMTK tools initialized.")

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
