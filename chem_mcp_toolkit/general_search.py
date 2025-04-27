import os
import logging
from tavily import TavilyClient

from utils.base_tool import BaseTool
from utils.errors import ChemMTKInputError, ChemMTKSearchFailError, ChemMTKToolProcessError, ChemMTKApiNotFoundError

from mcp_app import mcp


logger = logging.getLogger(__name__)


class WebSearch(BaseTool):
    name = "WebSearch"
    func_name = 'search_web'
    description = "Search the web for any questions and knowledge (including both general ones and domain-specific ones) and obtain concise summaries of the most relevant content. Input a specific question, returns a summary of the relevant content that answers the question."
    func_doc = ("question: str", "str")
    func_description = description
    examples = [
        {'input': 'What is the boiling point of water?', 'output': 'The boiling point of water at sea level is 100°C (212°F).'},
    ]

    def __init__(self, tavily_api_key: str, init=True, interface='text'):
        assert tavily_api_key is not None
        self.client = TavilyClient(api_key=tavily_api_key)
        super().__init__(init, interface=interface)

    def _run_base(self, query: str, *args, **kwargs) -> str:
        response = self.client.search(query, search_depth='advanced', include_answer=True)
        answer = response['answer']
        return answer
    

search_web = None

@mcp.tool()
async def search_web(query: str):
    """Search the web for any questions and knowledge (including both general ones and domain-specific ones) and obtain concise summaries of the most relevant content.
    
    Args:
        query: The search query.
    Returns:
        str: The summaries of related content.
    """
    tavily_api_key = os.getenv("TAVILY_API_KEY", None)
    if tavily_api_key is None:
        raise ChemMTKApiNotFoundError("Cannot find the API key for Tavily. Please set the TAVILY_API_KEY environment variable.")
    try:
        if search_web is None:
            search_web = WebSearch(tavily_api_key=tavily_api_key)
        return await search_web.run(query)
    except Exception as e:
        raise ChemMTKSearchFailError(f"Error searching the web: {e}") from e


wikipedia_wrapper = None


@mcp.tool()
async def search_wikipedia(query: str):
    """Search Wikipedia for a given query.
    
    Args:
        query: The search query.
    Returns:
        str: The summaries of related content.
    """
    try:
        if wikipedia_wrapper is None:
            import wikipedia
            global wikipedia_wrapper
            wikipedia_wrapper = wikipedia

        search_results = wikipedia_wrapper.search(query)
        if search_results:
            page = wikipedia_wrapper.page(search_results[0])
            return page.url
        else:
            raise ChemMTKSearchFailError("No search results found.")
    except Exception as e:
        raise ChemMTKSearchFailError(f"Error searching Wikipedia: {e}")

    


if __name__ == "__main__":
    # Initialize and run the server
    mcp.run(transport='stdio')
