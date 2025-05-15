---
categories:
- Molecule
description: 'Predict whether a molecule can cause 20 different side effects, along
  with the probabilities of each side effect. The side effects are: (1) Blood and
  lymphatic system disorders; (2) Cardiac disorders; (3) Congenital, familial and
  genetic disorders; (4) Ear and labyrinth disorders; (5) Endocrine disorders; (6)
  Eye disorders; (7) Gastrointestinal disorders; (8) Hepatobiliary disorders; (9)
  Immune system disorders; (10) Metabolism and nutrition disorders; (11) Musculoskeletal
  and connective tissue disorders; (12) Neoplasms benign, malignant and unspecified
  (incl cysts and polyps); (13) Nervous system disorders; (14) Pregnancy, puerperium
  and perinatal conditions; (15) Psychiatric disorders; (16) Renal and urinary disorders;
  (17) Reproductive system and breast disorders; (18) Respiratory, thoracic and mediastinal
  disorders; (19) Skin and subcutaneous tissue disorders; (20) Vascular disorders.'
draft: false
tags:
- Molecular Properties
- Neural Networks
title: SideEffectPredictor (predict_side_effect)
weight: 2

---
<div style="display: flex; flex-wrap: wrap; gap: 0.75rem; align-items: center;">
  {{< badge >}}Version: 0.1.0{{< /badge >}}
  {{< badge >}}Last Update: Unknown{{< /badge >}}
  {{< badge >}}MCP Support{{< /badge >}}
  {{< badge >}}Python Calling Support{{< /badge >}}
</div>
{{< lead >}}
**Predict whether a molecule can cause 20 different side effects, along with the probabilities of each side effect. The side effects are: (1) Blood and lymphatic system disorders; (2) Cardiac disorders; (3) Congenital, familial and genetic disorders; (4) Ear and labyrinth disorders; (5) Endocrine disorders; (6) Eye disorders; (7) Gastrointestinal disorders; (8) Hepatobiliary disorders; (9) Immune system disorders; (10) Metabolism and nutrition disorders; (11) Musculoskeletal and connective tissue disorders; (12) Neoplasms benign, malignant and unspecified (incl cysts and polyps); (13) Nervous system disorders; (14) Pregnancy, puerperium and perinatal conditions; (15) Psychiatric disorders; (16) Renal and urinary disorders; (17) Reproductive system and breast disorders; (18) Respiratory, thoracic and mediastinal disorders; (19) Skin and subcutaneous tissue disorders; (20) Vascular disorders.**
{{< /lead >}}

**Example**

Input:
```yaml
smiles: 'CC1=CC(C)=C(NC(=O)CN(CC(=O)O)CC(=O)O)C(C)=C1Br'
```

Text Input (used for the `run_text` function in the Python calling mode):
```yaml
smiles: 'CC1=CC(C)=C(NC(=O)CN(CC(=O)O)CC(=O)O)C(C)=C1Br'
```

Output:
```yaml
side_effect: "The probabilities of the compound to cause different side effects are as follows:Blood and lymphatic system disorders: 11.29%, which means it's unlikely to cause the side effect.\nCardiac disorders: 10.92%, which means it's unlikely to cause the side effect.\nCongenital, familial and genetic disorders: 11.98%, which means it's unlikely to cause the side effect.\nEar and labyrinth disorders: 8.48%, which means it's unlikely to cause the side effect.\nEndocrine disorders: 4.16%, which means it's unlikely to cause the side effect.\nEye disorders: 15.19%, which means it's unlikely to cause the side effect.\nGastrointestinal disorders: 57.00%, which means it's likely to cause the side effect.\nHepatobiliary disorders: 9.62%, which means it's unlikely to cause the side effect.\nImmune system disorders: 10.14%, which means it's unlikely to cause the side effect.\nMetabolism and nutrition disorders: 15.41%, which means it's unlikely to cause the side effect.\nMusculoskeletal and connective tissue disorders: 10.77%, which means it's unlikely to cause the side effect.\nNeoplasms benign, malignant and unspecified (incl cysts and polyps): 4.92%, which means it's unlikely to cause the side effect.\nNervous system disorders: 34.37%, which means it's unlikely to cause the side effect.\nPregnancy, puerperium and perinatal conditions: 3.32%, which means it's unlikely to cause the side effect.\nPsychiatric disorders: 8.06%, which means it's unlikely to cause the side effect.\nRenal and urinary disorders: 10.64%, which means it's unlikely to cause the side effect.\nReproductive system and breast disorders: 4.59%, which means it's unlikely to cause the side effect.\nRespiratory, thoracic and mediastinal disorders: 16.48%, which means it's unlikely to cause the side effect.\nSkin and subcutaneous tissue disorders: 53.97%, which means it's likely to cause the side effect.\nVascular disorders: 18.45%, which means it's unlikely to cause the side effect.\nNote that the results are predicted by a neural network model and may not be accurate. You may use other tools or resources to obtain more reliable results if needed."
```

## Usage

The tool supports both [MCP mode](#mcp-mode) and [Python calling mode](#python-calling-mode).



### MCP Mode

Configure your MCP client following its instructions with something like:
```JSON
{
    "command": "/ABSTRACT/PATH/TO/uv",  // Use `which uv` to get its path
    "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run", "-m", "chemmcp.tools.side_effect_predictor"],
    "toolCallTimeoutMillis": 300000,
    "env": {}
}
```

### Python Calling Mode

```python
import os
from chemmcp.tools import SideEffectPredictor

# Initialize the tool
tool = SideEffectPredictor()

# The tool has two alternative ways to run:
# 1. Run with separate input domains (recommended)
output = tool.run_code(
    smiles='CC1=CC(C)=C(NC(=O)CN(CC(=O)O)CC(=O)O)C(C)=C1Br'
)
# 2. Run with text-only input
output = tool.run_text(
    smiles='CC1=CC(C)=C(NC(=O)CN(CC(=O)O)CC(=O)O)C(C)=C1Br'
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

### Text Input
Used in the `run_text` function in the Python calling mode.
| Name | Type | Default | Description |
| --- | --- | --- | --- |
| smiles | str | N/A | SMILES string of the molecule. |

### Output
The output is the same in both input cases.
| Name | Type | Description |
| --- | --- | --- |
| side_effect | str | The probabilities of the compound to cause different side effects. |

### Envs
No required environment variables for this tool.
