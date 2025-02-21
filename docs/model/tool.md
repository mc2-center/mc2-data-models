A **Tool** entry documents software applications, utilities, or platforms designed to perform specific operations in research or technical workflows. These tools can include anything from data analysis programs to machine learning frameworks and bioinformatics utilities. Accurate documentation ensures that tools are easily discoverable, traceable, and appropriately integrated into research processes.

### **Why You Should Contribute Tool Entries**
Contributing tool entries supports collaboration, reproducibility, and innovation by helping researchers and teams discover and use essential software for data analysis, visualization, and processing. Accurate documentation ensures tools are accessible and impactful across projects.

#### **Who Should Be Contributing Tool Entries?**
1. **Developers** – Share tools to boost adoption and integration.
2. **Researchers** – Document tools for replicability and collaboration.
3. **Data Scientists** – Ensure tools critical to workflows are recorded.
4. **Project Managers** – Maintain visibility of tools used across teams.
5. **IT and Support Specialists** – Provide operational details and support resources.

## Download Template
You can download the [ToolView CSV template](https://github.com/mc2-center/data-models/raw/main/templates/ToolView.csv) to streamline data entry.


## Example Data Entry
The table below includes sample values to demonstrate proper attribute usage.

| **Attribute**                 | **Example Value**                                                                                                    |
|-------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| Tool Name                     | Adobe Photoshop                                                                                                         |
| Tool Description              | "A photo editing tool with advanced image manipulation features, supporting multiple file formats and layers."          |
| Tool Release Date             | January 5, 2022                                                                                                        |
| Tool Entity Role              | Developer, Maintainer                                                                                                  |
| Tool Entity Type              | Organization                                                                                                           |
| Tool Input Data               | DNA Sequence                                                                                                           |
| Tool Output Data              | Gene ID (NCBI): NM_001282392.1                                                                                         |
| Tool Grant Number             | CA209975                                                                                                                |
| Tool Documentation Url        | https://docs.example.com/tool-guide.html                                                                               |
| Tool Operating System         | Windows, MacOS                                                                                                          |
| Tool Version                  | 3.2.1                                                                                                                  |
| Tool View                     | Detail View                                                                                                             |
| Tool Pubmed Id                | 26760201                                                                                                                |
| Tool License                  | Apache-2.0                                                                                                              |

## Full Field Reference

Below is the full field reference table with attributes and their descriptions.

[⤓ Download template](https://github.com/mc2-center/data-models/raw/main/templates/ToolView.csv)

{{ read_csv('tool/template.csv') }}
