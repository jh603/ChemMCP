---
categories:
- Molecule
description: Canonicalize a molecular SMILES string.
draft: false
tags:
- SMILES
- Canonicalization
- RDKit
title: SmilesCanonicalization (canonicalize_smiles)
weight: 2

---
<div style="display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;">
  {{< badge >}}Version: 0.1.0{{< /badge >}}
  {{< badge >}}Last Update: 2025/05/17{{< /badge >}}
  {{< badge >}}MCP Support{{< /badge >}}
  {{< badge >}}Python Calling Support{{< /badge >}}
</div>
{{< lead >}}
**Canonicalize a molecular SMILES string.**
{{< /lead >}}

**Example**

Input:
```yaml
smiles: 'C(O)C'
isomeric: True
kekulization: True
keep_atom_map: False
```

Text Input (used for the `run_text` function in the Python calling mode):
```yaml
smiles: 'C(O)C'
```

Output:
```yaml
canonical_smiles: 'CCO'
```

## Usage

The tool supports both [MCP mode](#mcp-mode) and [Python calling mode](#python-calling-mode).



### MCP Mode

Configure your MCP client following its instructions with something like:
```JSON
{
    "command": "/ABSTRACT/PATH/TO/uv",  // Use `which uv` to get its path
    "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run", "--tools", "SmilesCanonicalization"],
    "toolCallTimeoutMillis": 300000,
    "env": {}
}
```

### Python Calling Mode

```python
import os
from chemmcp.tools import SmilesCanonicalization

# Initialize the tool
tool = SmilesCanonicalization()

# The tool has two alternative ways to run:
# 1. Run with separate input domains (recommended)
output = tool.run_code(
    smiles='C(O)C'
    isomeric=True
    kekulization=True
    keep_atom_map=False
)
# 2. Run with text-only input
output = tool.run_text(
    smiles='C(O)C'
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
| smiles | str | N/A | SMILES string of the molecule. |
| isomeric | bool | N/A | Whether to include isomeric information. Default is True. |
| kekulization | bool | N/A | Whether to perform kekulization. Default is True. |
| keep_atom_map | bool | N/A | Whether to keep atom mapping numbers, if any. Default is True. |

### Text Input
Used in the `run_text` function in the Python calling mode.
| Name | Type | Default | Description |
| --- | --- | --- | --- |
| smiles | str | N/A | SMILES string of the molecule. |

### Output
The output is the same in both input cases.
| Name | Type | Description |
| --- | --- | --- |
| canonical_smiles | str | Canonicalized SMILES string. |

### Envs
No required environment variables for this tool.
