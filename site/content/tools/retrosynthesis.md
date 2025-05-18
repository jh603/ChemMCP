---
categories:
- Reaction
description: Conduct single-step retrosynthesis. Given the product(s), predict multiple
  sets of potential reactants, along with their confidence.
draft: false
tags:
- Retrosynthesis
- Reaction Prediction
title: Retrosynthesis (do_retrosynthesis)
weight: 2

---
<div style="display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;">
  {{< badge >}}Version: 0.1.0{{< /badge >}}
  {{< badge >}}Last Update: 2025/05/17{{< /badge >}}
  {{< badge >}}MCP Support{{< /badge >}}
  {{< badge >}}Python Calling Support{{< /badge >}}
</div>
{{< lead >}}
**Conduct single-step retrosynthesis. Given the product(s), predict multiple sets of potential reactants, along with their confidence.**
{{< /lead >}}

**Example**

Input:
```yaml
product_smiles: 'CCO'
```

Text Input (used for the `run_text` function in the Python calling mode):
```yaml
product_smiles: 'CCO'
```

Output:
```yaml
reactants_and_confidence: 'There are 13 possible sets of reactants for the given product:\n1.\tReactants: C1CCOC1.CCNC(=O)c1cccn1C.[Li][AlH4]\tConfidence: 1.0\n2.\tReactants: CCN.CCO.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n3.\tReactants: CCN.CO.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n4.\tReactants: CCN.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n5.\tReactants: CCN.CCO.Cn1cccc1C=O.O.[BH4-].[Na+]\tConfidence: 1.0\n6.\tReactants: CCN.CO.Cn1cccc1C=O.O.[BH4-].[Na+]\tConfidence: 1.0\n7.\tReactants: C1CCOC1.CCN.Cn1cccc1C=O.[BH4-].[Na+]\tConfidence: 1.0\n8.\tReactants: CCN.Cl.Cn1cccc1C=O\tConfidence: 0.938\n9.\tReactants: CCN.Cn1cccc1C=O\tConfidence: 0.917\n10.\tReactants: CCN.Cl.Cn1cccc1C=O\tConfidence: 0.841\n11.\tReactants: C1CCOC1.CCN.Cn1cccc1C=O\tConfidence: 0.797\n12.\tReactants: C1CCOC1.CCN.CO.Cn1cccc1C=O\tConfidence: 0.647\n13.\tReactants: C1CCOC1.CC(=O)NCc1cccn1C.[Li][AlH4]\tConfidence: 1.0\n'
```

## Usage

The tool supports both [MCP mode](#mcp-mode) and [Python calling mode](#python-calling-mode).

### Environment Variables
This tool requires the following environment variables:
- **RXN4CHEM_API_KEY**: The API key for IBM RXN4Chem.


### MCP Mode

Configure your MCP client following its instructions with something like:
```JSON
{
    "command": "/ABSTRACT/PATH/TO/uv",  // Use `which uv` to get its path
    "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run", "--tools", "Retrosynthesis"],
    "toolCallTimeoutMillis": 300000,
    "env": {
        "RXN4CHEM_API_KEY": "VALUE_TO_BE_SET"
    }
}
```

### Python Calling Mode

```python
import os
from chemmcp.tools import Retrosynthesis

# Set the environment variables
os.environ['RXN4CHEM_API_KEY'] = 'VALUE_TO_BE_SET'

# Initialize the tool
tool = Retrosynthesis()

# The tool has two alternative ways to run:
# 1. Run with separate input domains (recommended)
output = tool.run_code(
    product_smiles='CCO'
)
# 2. Run with text-only input
output = tool.run_text(
    product_smiles='CCO'
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
| product_smiles | str | N/A | The SMILES of the product. |

### Text Input
Used in the `run_text` function in the Python calling mode.
| Name | Type | Default | Description |
| --- | --- | --- | --- |
| product_smiles | str | N/A | The SMILES of the product. |

### Output
The output is the same in both input cases.
| Name | Type | Description |
| --- | --- | --- |
| reactants_and_confidence | str | The SMILES of the reactants and the confidence. |

### Envs
| Name | Description |
| --- | --- |
| RXN4CHEM_API_KEY | The API key for IBM RXN4Chem. |
