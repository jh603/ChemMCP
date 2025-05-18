import json
from typing import List


def generate_tools_env_json(tools: List[str], env_path: str='site/data/tool_envs/all_tool_envs.json'):
    with open(env_path, 'r') as f:
        tools_envs_dict = json.load(f)

    final_tools_envs_dict = {}
    for tool_name in tools:
        tool_envs_dict = tools_envs_dict[tool_name]
        final_tools_envs_dict.update(tool_envs_dict)
            
    if '__llms__' in tools_envs_dict:
        tools_envs_dict.pop('__llms__')
        tools_envs_dict['LLM_MODEL_NAME'] = 'The name of the LLM to use. See [LiteLLM](https://docs.litellm.ai/docs/#basic-usage) for more details.'
        tools_envs_dict['OTHER_LLM_CREDENTIALS'] = 'Other LLM credentials are required to be set in the `env` field. See [LiteLLM](https://docs.litellm.ai/docs/#basic-usage) for more details.'

    all_tools = list(tools_envs_dict.keys())

    return final_tools_envs_dict, all_tools


def generate_markdown_template(tools: List[str]):

    envs, all_tools = generate_tools_env_json(tools, env_path='site/data/tool_envs/all_tool_envs.json')

    if len(tools) == len(all_tools):
        select_tool_command = []
    else:
        select_tool_command = ["--tools"] + tools

    config_template = {
        "command": "/ABSTRACT/PATH/TO/uv",
        "args": ["--directory", "/ABSTRACT/PATH/TO/ChemMCP", "run"] + select_tool_command,
        "toolCallTimeoutMillis": 300000,
        "env": envs
    }

    config_str = json.dumps(config_template, indent=4)

    markdown = f"""Here is the configuration for using your chosen {len()} tools in ChemMCP:

```JSON
{config_str}
```
"""
    
    return markdown
