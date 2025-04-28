from . import molecule_ops
from . import molecule_info
from . import name_conversion
from . import molecule_description
from . import property_prediction
from . import pubchem_search

from .mcp_app import mcp


# if __name__ == "__main__":
#     # Initialize and run the server
#     mcp.run(transport='stdio')


# build a Starlette/uvicorn app
app = mcp.sse_app()

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8001, log_level="info")
