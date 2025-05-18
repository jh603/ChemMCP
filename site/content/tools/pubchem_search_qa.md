---
categories:
- Molecule
description: Answer questions about molecules/compounds based on the information from
  PubChem, one of the most comprehensive database of chemical molecules and their
  activities.
draft: false
tags:
- PubChem
- Molecule Information
- Molecular Properties
title: PubchemSearchQA (search_pubchem_qa)
weight: 2

---
<div style="display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;">
  {{< badge >}}Version: 0.1.0{{< /badge >}}
  {{< badge >}}Last Update: 2025/05/17{{< /badge >}}
  {{< badge >}}MCP Support{{< /badge >}}
  {{< badge >}}Python Calling Support{{< /badge >}}
</div>
{{< lead >}}
**Answer questions about molecules/compounds based on the information from PubChem, one of the most comprehensive database of chemical molecules and their activities.**
{{< /lead >}}

**Example 3**

Input:
```yaml
representation_name: 'Name'
representation: 'alcohol'
question: 'What properties do this molecule have?'
```

Text Input (used for the `run_text` function in the Python calling mode):
```yaml
representation_name_and_representation_and_question: 'Name: alcohol   Questions: What properties do this molecule have?'
```

Output:
```yaml
answer: 'This molecule has the following properties: [...]'
```

## Usage

The tool supports both [MCP mode](#mcp-mode) and [Python calling mode](#python-calling-mode).

### Environment Variables
This tool requires the following environment variables:
- **LLM_MODEL_NAME**: The name of the LLM to use. See [LiteLLM](https://github.com/Lightning-AI/litellm) for more details.
- Other LLM credentials are required to be set in the `env` field. See [LiteLLM](https://github.com/Lightning-AI/litellm) for more details.


### MCP Mode

Configure your MCP client following its instructions with something like:
```JSON
{
    "command": "/ABSTRACT/PATH/TO/uv",  // Use `which uv` to get its path
    "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run", "--tools", "PubchemSearchQA"],
    "toolCallTimeoutMillis": 300000,
    "env": {
        "LLM_MODEL_NAME": "VALUE_TO_BE_SET"
        // Add required LLM credentials
        // ...
    }
}
```

### Python Calling Mode

```python
import os
from chemmcp.tools import PubchemSearchQA

# Set the environment variables
os.environ['LLM_MODEL_NAME'] = 'VALUE_TO_BE_SET'
# Also add LLM credentials required by the LLM model, such as `OPENAI_API_KEY`

# Initialize the tool
tool = PubchemSearchQA()

# The tool has two alternative ways to run:
# 1. Run with separate input domains (recommended)
output = tool.run_code(
    representation_name='SMILES'
    representation='CCO'
    question='What properties do this molecule have?'
)
# 2. Run with text-only input
output = tool.run_text(
    representation_name_and_representation_and_question='SMILES: CCO   Questions: What properties do this molecule have?'
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
| representation_name | str | N/A | The representation name, can be "smiles", "iupac", or "name" (chemical's common name). |
| representation | str | N/A | The representation of the molecule/compound, corresponding to the representation_name used. |
| question | str | N/A | The question about the molecule/compound. |

### Text Input
Used in the `run_text` function in the Python calling mode.
| Name | Type | Default | Description |
| --- | --- | --- | --- |
| representation_name_and_representation_and_question | str | N/A | The representation name and representation of the molecule/compound, e.g., "SMILES: <SMILES>", "IUPAC: <IUPAC name>", or "Name: <common name>". Followed by "Question: <your question about the molecule/compound>". |

### Output
The output is the same in both input cases.
| Name | Type | Description |
| --- | --- | --- |
| answer | str | The answer to the question based on the PubChem page. |

### Envs
| Name | Description |
| --- | --- |
| LLM_MODEL_NAME | The name of the LLM to use. See [LiteLLM](https://github.com/Lightning-AI/litellm) for more details. |
| LLM credentials | Other LLM credentials are required to be set in the `env` field. See [LiteLLM](https://github.com/Lightning-AI/litellm) for more details. |
