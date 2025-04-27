import molecule_ops
import molecule_info
import name_conversion

from mcp_app import mcp



if __name__ == "__main__":
    # Initialize and run the server
    mcp.run(transport='stdio')
