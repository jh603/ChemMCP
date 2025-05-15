---
categories:
- General
description: Search the web for any questions and knowledge and obtain a concise answer
  based on thesearch results.
draft: false
tags:
- Web
- LLMs
- Neural Networks
title: WebSearch (search_web)
weight: 2

---
<div style="display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;">
  {{< badge >}}Version: 0.1.0{{< /badge >}}
  {{< badge >}}Last Update: Unknown{{< /badge >}}
  {{< badge >}}MCP Support{{< /badge >}}
  {{< badge >}}Python Calling Support{{< /badge >}}
</div>
{{< lead >}}
**Search the web for any questions and knowledge and obtain a concise answer based on thesearch results.**
{{< /lead >}}

**Example**

Input:
```yaml
query: 'What is the boiling point of water?'
```

Text Input (used for the `run_text` function in the Python calling mode):
```yaml
query: 'What is the boiling point of water?'
```

Output:
```yaml
result: 'The boiling point of water at sea level is 100°C (212°F).'
```

## Usage

The tool supports both [MCP mode](#mcp-mode) and [Python calling mode](#python-calling-mode).

### Environment Variables
This tool requires the following environment variables:
- **TAVILY_API_KEY**: The API key for [Tavily](https://tavily.com/).


### MCP Mode

Configure your MCP client following its instructions with something like:
```JSON
{
    "command": "/ABSTRACT/PATH/TO/uv",  // Use `which uv` to get its path
    "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run", "-m", "chemmcp.tools.web_search"],
    "toolCallTimeoutMillis": 300000,
    "env": {
        "TAVILY_API_KEY": "VALUE_TO_BE_SET"
    }
}
```

### Python Calling Mode

```python
import os
from chemmcp.tools import WebSearch

# Set the environment variables
os.environ['TAVILY_API_KEY'] = 'VALUE_TO_BE_SET'

# Initialize the tool
tool = WebSearch()

# The tool has two alternative ways to run:
# 1. Run with separate input domains (recommended)
output = tool.run_code(
    query='What is the boiling point of water?'
)
# 2. Run with text-only input
output = tool.run_text(
    query='What is the boiling point of water?'
)
```


Each tool in ChemMCP has two ways to run:
- **`run_code`** (recommended): The inputs contain one or more domains, each of which can be a str, an int, a float, etc.
- **`run_text`**: The inputs are a single string in a specific format. The tool will parse the string to extract the input domains. This is useful in scenarios where an agent framework calls tools only with text input.
The output is the same in both cases.

For the input and output domains, please refer to the tool's [signature](#tool-signature).

## Tool Signature



### Input
Used in the MCP mode, as well as the `run_code` function in the Python calling mode.
| Name | Type | Default | Description |
| --- | --- | --- | --- |
| query | str | N/A | The search query. |

### Text Input
Used in the `run_text` function in the Python calling mode.
| Name | Type | Default | Description |
| --- | --- | --- | --- |
| query | str | N/A | The search query. |

### Output
The output is the same in both input cases.
| Name | Type | Description |
| --- | --- | --- |
| result | str | The answer to the search query summarized by Tavily's LLM. |

### Envs
| Name | Description |
| --- | --- |
| TAVILY_API_KEY | The API key for [Tavily](https://tavily.com/). |
