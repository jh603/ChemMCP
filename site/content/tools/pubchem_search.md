---
categories:
- Molecule
description: Search for molecule/compound information on PubChem, one of the most
  comprehensive database of chemical molecules and their activities.
draft: false
tags:
- PubChem
- Molecule Information
- Molecular Properties
title: PubchemSearch (search_pubchem)
weight: 2

---
<div style="display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;">
  {{< badge >}}Version: 0.1.0{{< /badge >}}
  {{< badge >}}Last Update: Unknown{{< /badge >}}
  {{< badge >}}No MCP Support{{< /badge >}}
  {{< badge >}}Python Calling Support{{< /badge >}}
</div>
{{< lead >}}
**Search for molecule/compound information on PubChem, one of the most comprehensive database of chemical molecules and their activities.**
{{< /lead >}}

**Example 3**

Input:
```yaml
representation_name: 'Name'
representation: 'alcohol'
```

Text Input (used for the `run_text` function in the Python calling mode):
```yaml
representation_name_and_representation: 'Name: alcohol'
```

Output:
```yaml
compound_doc: '# 1 Names and Identifiers\nSection Description: Chemical names, synonyms, identifiers, and descriptors.\n\n## 1.1 Record Description\nSection Description: Summary Information\n\nEthanol with a small amount of an adulterant added so as to be unfit for use as a beverage. [...]'
```

## Usage

The tool supports [Python calling mode](#python-calling-mode).



### MCP Mode

This tool does not support MCP mode.

### Python Calling Mode

```python
import os
from chemmcp.tools import PubchemSearch

# Initialize the tool
tool = PubchemSearch()

# The tool has two alternative ways to run:
# 1. Run with separate input domains (recommended)
output = tool.run_code(
    representation_name='SMILES'
    representation='CCO'
)
# 2. Run with text-only input
output = tool.run_text(
    representation_name_and_representation='SMILES: CCO'
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
| representation_name | str | N/A | The representation name, can be "SMILES", "IUPAC", or "Name" (chemical's common name). |
| representation | str | N/A | The representation of the molecule/compound, corresponding to the representation_name used. |

### Text Input
Used in the `run_text` function in the Python calling mode.
| Name | Type | Default | Description |
| --- | --- | --- | --- |
| representation_name_and_representation | str | N/A | The representation name and representation of the molecule/compound, e.g., "SMILES: <SMILES>", "IUPAC: <IUPAC name>", or "Name: <common name>". |

### Output
The output is the same in both input cases.
| Name | Type | Description |
| --- | --- | --- |
| compound_doc | str | The document of the molecule/compound in a markdown format. |

### Envs
No required environment variables for this tool.
