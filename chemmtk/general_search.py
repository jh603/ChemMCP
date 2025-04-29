import os
import logging
import argparse

from tavily import TavilyClient

from .utils.base_tool import BaseTool, register_mcp_tool
from .utils.errors import ChemMTKInputError, ChemMTKSearchFailError, ChemMTKToolProcessError, ChemMTKApiNotFoundError
from .mcp_app import mcp


logger = logging.getLogger(__name__)


@register_mcp_tool(mcp)
class WebSearch(BaseTool):
    name = "WebSearch"
    func_name = 'search_web'
    description = "Search the web for any questions and knowledge (including both general ones and domain-specific ones) and obtain a concise answer of the search results."
    text_input = [("query", "str", "The search query.")]
    code_input = [("query", "str", "The search query.")]
    tool_output = [("result", "str", "The summaries of related content.")]
    examples = [
        {'input': 'What is the boiling point of water?', 'output': 'The boiling point of water at sea level is 100°C (212°F).'},
    ]

    def __init__(self, tavily_api_key: str = None, init=True, interface='text'):
        if tavily_api_key is None:
            tavily_api_key = os.getenv("TAVILY_API_KEY", None)
        if tavily_api_key is None:
            raise ChemMTKApiNotFoundError("Cannot find the API key for Tavily. Please set the TAVILY_API_KEY environment variable.")
        self.client = TavilyClient(api_key=tavily_api_key)
        super().__init__(init, interface=interface)

    def _run_base(self, query: str) -> str:
        logger.info("Running WebSearch with query: %s", query)
        response = self.client.search(query, search_depth='advanced', include_answer=True)
        try:
            answer = response['answer']
        except KeyError as e:
            raise ChemMTKSearchFailError(f"Error searching the web: {e}") from e
        return answer


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the MCP server.")
    parser.add_argument('--sse', action='store_true', help="Run the server with SSE (Server-Sent Events) support.")
    args = parser.parse_args()

    if args.sse:
        # build a Starlette/uvicorn app
        app = mcp.sse_app()
        import uvicorn
        uvicorn.run(app, host="127.0.0.1", port=8001)
    else:
        # Run the MCP server with standard input/output
        mcp.run(transport='stdio')
