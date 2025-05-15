---
categories:
- Molecule
description: Generate a textual description of the molecule from its SMILES representation
  with MolT5. This tool uses neural networks to generate descriptions, which may not
  be accurate or correct. Please first try other tools that provide accurate and authoritative
  information, and only use this one as the last resort.
draft: false
tags:
- Molecular Description
- Neural Networks
title: MoleculeCaptioner (generate_molecule_caption)
weight: 2

---
<div style="display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;">
  {{< badge >}}Version: 0.1.0{{< /badge >}}
  {{< badge >}}Last Update: 2025/05/14{{< /badge >}}
  {{< badge >}}MCP Support{{< /badge >}}
  {{< badge >}}Python Calling Support{{< /badge >}}
</div>
{{< lead >}}
**Generate a textual description of the molecule from its SMILES representation with MolT5. This tool uses neural networks to generate descriptions, which may not be accurate or correct. Please first try other tools that provide accurate and authoritative information, and only use this one as the last resort.**
{{< /lead >}}

**Example**

Input:
```yaml
smiles: 'CCO'
```

Text Input (used for the `run_text` function in the Python calling mode):
```yaml
smiles: 'CCO'
```

Output:
```yaml
description: 'The molecule is an ether in which the oxygen atom is linked to two ethyl groups. It has a role as an inhalation anaesthetic, a non-polar solvent and a refrigerant. It is a volatile organic compound and an ether.\n\nNote: This is a generated description and may not be accurate. Please double check the result.'
```

## Usage

The tool supports both [MCP mode](#mcp-mode) and [Python calling mode](#python-calling-mode).



### MCP Mode

Configure your MCP client following its instructions with something like:
```JSON
{
    "command": "/ABSTRACT/PATH/TO/uv",  // Use `which uv` to get its path
    "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run", "-m", "chemmcp.tools.molecule_captioner"],
    "toolCallTimeoutMillis": 300000,
    "env": {}
}
```

### Python Calling Mode

```python
import os
from chemmcp.tools import MoleculeCaptioner

# Initialize the tool
tool = MoleculeCaptioner()

# The tool has two alternative ways to run:
# 1. Run with separate input domains (recommended)
output = tool.run_code(
    smiles='CCO'
)
# 2. Run with text-only input
output = tool.run_text(
    smiles='CCO'
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
| smiles | str | N/A | SMILES representation of the molecule. |

### Text Input
Used in the `run_text` function in the Python calling mode.
| Name | Type | Default | Description |
| --- | --- | --- | --- |
| smiles | str | N/A | SMILES representation of the molecule. |

### Output
The output is the same in both input cases.
| Name | Type | Description |
| --- | --- | --- |
| description | str | Textual description of the molecule. |

### Envs
No required environment variables for this tool.
