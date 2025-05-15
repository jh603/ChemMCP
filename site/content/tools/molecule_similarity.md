---
categories:
- Molecule
description: Get the Tanimoto similarity of two molecules. It Can also be used to
  check if two molecules are identical.
draft: false
tags:
- Molecular Property
- RDKit
title: MoleculeSimilarity (cal_molecule_similarity)
weight: 2

---
<div style="display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;">
  {{< badge >}}Version: 0.1.0{{< /badge >}}
  {{< badge >}}Last Update: 2025/05/14{{< /badge >}}
  {{< badge >}}MCP Support{{< /badge >}}
  {{< badge >}}Python Calling Support{{< /badge >}}
</div>
{{< lead >}}
**Get the Tanimoto similarity of two molecules. It Can also be used to check if two molecules are identical.**
{{< /lead >}}

**Example 2**

Input:
```yaml
smiles1: 'CCO'
smiles2: 'C(O)C'
```

Text Input (used for the `run_text` function in the Python calling mode):
```yaml
smiles_pair: 'CCO;C(O)C'
```

Output:
```yaml
similarity: 'Input Molecules Are Identical'
```

## Usage

The tool supports both [MCP mode](#mcp-mode) and [Python calling mode](#python-calling-mode).



### MCP Mode

Configure your MCP client following its instructions with something like:
```JSON
{
    "command": "/ABSTRACT/PATH/TO/uv",  // Use `which uv` to get its path
    "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run", "-m", "chemmcp.tools.molecule_similarity"],
    "toolCallTimeoutMillis": 300000,
    "env": {}
}
```

### Python Calling Mode

```python
import os
from chemmcp.tools import MoleculeSimilarity

# Initialize the tool
tool = MoleculeSimilarity()

# The tool has two alternative ways to run:
# 1. Run with separate input domains (recommended)
output = tool.run_code(
    smiles1='CCO'
    smiles2='CCN'
)
# 2. Run with text-only input
output = tool.run_text(
    smiles_pair='CCO;CCN'
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
| smiles1 | str | N/A | SMILES string of the first molecule |
| smiles2 | str | N/A | SMILES string of the second molecule. |

### Text Input
Used in the `run_text` function in the Python calling mode.
| Name | Type | Default | Description |
| --- | --- | --- | --- |
| smiles_pair | str | N/A | SMILES strings of the two molecules, separated by a semicolon. |

### Output
The output is the same in both input cases.
| Name | Type | Description |
| --- | --- | --- |
| similarity | str | Tanimoto similarity score and similarity description |

### Envs
No required environment variables for this tool.
