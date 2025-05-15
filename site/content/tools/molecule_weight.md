---
categories:
- Molecule
description: Calculate molecular weight.
draft: false
tags:
- Molecular Property
- RDKit
title: MoleculeWeight (cal_molecular_weight)
weight: 2

---
<div style="display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;">
  {{< badge >}}Version: 0.1.0{{< /badge >}}
  {{< badge >}}Last Update: Unknown{{< /badge >}}
  {{< badge >}}MCP Support{{< /badge >}}
  {{< badge >}}Python Calling Support{{< /badge >}}
</div>
{{< lead >}}
**Calculate molecular weight.**
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
weight: 46.041864812
```

## Usage

The tool supports both [MCP mode](#mcp-mode) and [Python calling mode](#python-calling-mode).



### MCP Mode

Configure your MCP client following its instructions with something like:
```JSON
{
    "command": "/ABSTRACT/PATH/TO/uv",  // Use `which uv` to get its path
    "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run", "-m", "chemmcp.tools.molecule_weight"],
    "toolCallTimeoutMillis": 300000,
    "env": {}
}
```

### Python Calling Mode

```python
import os
from chemmcp.tools import MoleculeWeight

# Initialize the tool
tool = MoleculeWeight()

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
| smiles | str | N/A | SMILES string of the molecule |

### Text Input
Used in the `run_text` function in the Python calling mode.
| Name | Type | Default | Description |
| --- | --- | --- | --- |
| smiles | str | N/A | SMILES string of the molecule |

### Output
The output is the same in both input cases.
| Name | Type | Description |
| --- | --- | --- |
| weight | float | Molecular weight of the molecule |

### Envs
No required environment variables for this tool.
