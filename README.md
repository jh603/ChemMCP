# ChemMCP

<img src="site/assets/img/icon_text_logo.png" alt="ChemMCP Banner" style="zoom:20%;" />

See our website and document at [osu-nlp-group.github.io/ChemMCP](https://osu-nlp-group.github.io/ChemMCP/).

## What is this?

ChemMCP is **a continuously updated collection of chemistry tools for LLMs and AI assistants**, compatible with the [Model Context Protocol (MCP)](https://modelcontextprotocol.org/). By integrating **powerful chemistry tools**, ChemMCP **transforms general AI models into chemistry coscientists** capable of performing complex molecular analysis, property prediction, and reaction synthesis tasks without requiring domain-specific training.

[MCP (Model Context Protocol)](https://modelcontextprotocol.io/introduction) is a framework that allows AI models to access external tools and resources through a standardized interface. ChemMCP leverages this architecture to bridge the gap between general-purpose AI models and specialized chemistry tools, enabling seamless integration of chemistry expertise into AI workflows.

Specifically, ChemMCP provides:

- **Professional Chemistry Tools**: ChemMCP's curated collection of chemistry tools empowers you and your AI to predict molecules and reactions, analyze chemical data, explore scientific knowledge, and much more. Discover the full range in our [tool directory](https://osu-nlp-group.github.io/ChemMCP/tools).
- **Effortless Integration with MCP**: Seamlessly supercharge your LLM clients—like Claude, GPT, and others—using ChemMCP. Thanks to the MCP, you can add powerful chemistry tools to your workflow in minutes, simply by including a JSON configuration file.
- **Python-Ready for Your Projects**: Bring ChemMCP's power directly into your Python projects! Effortlessly integrate our tools for data processing, agent building, or your custom applications, perfect for researchers, developers, and innovators alike.

We will continue to add and maintain tools in ChemMCP. **You are more than welcome to contribute, by maining existing tools or adding new tools!**

## Get Started

[This document](https://osu-nlp-group.github.io/ChemMCP/get-started/) will help you quickly set up and start using ChemMCP.

## Tool List

Based on the functions, the tools can be divided into general tools, molecule tools, and reaction tools:

- **General Tools**: Provide broad information retrieval and web searching.
- **Molecule Tools**: Offer various analyses, predictions, and conversions related to chemical compounds and their properties.
- **Reaction Tools**: Predict products of chemical reactions and suggest potential reactants for synthesizing given products.

Check the updated full list of tools [here](https://osu-nlp-group.github.io/ChemMCP/tools/).

## Citation

If ChemMCP is valuable to your research or development, please kindly cite our work.

```
@misc{yu2025chemmcp,
  author       = {Botao Yu and Huan Sun},
  title        = {ChemMCP: A Chemistry MCP Toolkit},
  year         = {2025},
  url          = {https://github.com/OSU-NLP-Group/ChemMCP},
  note         = {2025-04-28-01},
}

@article{yu2024chemtoolagent,
    title={ChemToolAgent: The Impact of Tools on Language Agents for Chemistry Problem Solving},
    author={Botao Yu and Frazier N. Baker and Ziru Chen and Garrett Herb and Boyu Gou and Daniel Adu-Ampratwum and Xia Ning and Huan Sun},
    journal={arXiv preprint arXiv:2411.07228},
    year={2024}
}
```

