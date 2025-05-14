---
title: "WebSearch (search_web)"
categories: ["General"]
tags: ["Web", "LLMs", "Neural Networks"]
description: "Search the web for any questions and knowledge (including both general ones and domain-specific ones) and obtain a concise answer of the search results."
weight: 2
draft: false
---

Version: 0.1.0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Last Update: 2025/05/10

{{< lead >}}
**Search the web for any questions and knowledge (including both general ones and domain-specific ones) and obtain a concise answer of the search results.** 
{{< /lead >}}

**Example**

Input:
```yaml
query: "What is the boiling point of water?"
```

Output:
```
"The boiling point of water at sea level is 100°C (212°F)."
```

## Usage

This tool supports both [**MCP usage**](#mcp-usage) and [**Python calling**](#python-calling).

### MCP Usage

### Python Calling


## Tool Reference

### Input

| Name | Type | Default | Description |
| ----- | --- | --- | --- |
| query   | str  | N/A | The search query. |

### Output

| Name | Type | Default |Description |
| ----- | --- | --- | --- |
| result   | str | N/A | The answer to the search query summarized by Tavily's LLM. |

### Required Env

| Name  | Description |
| ----- | --- |
| TAVILY_API_KEY | The API key for [Tavily](https://tavily.com/).  |

