A **Dataset Sharing Plan (DataDSP)** is a structured framework used to document the sharing, storage, and usage details of datasets within  MC<sup>2</sup> Center-supported Synapse projects. These plans help ensure datasets are traceable, well-organized, and compliant with regulations regarding data sharing, accessibility, and ethical use. By defining attributes such as dataset names, assays, species, and sharing permissions, the model facilitates efficient data management and collaboration across research projects.

This page outlines the key attributes required to create a DataDSP, guiding users on how to structure and share datasets while adhering to best practices. It also demonstrates how attributes like planned upload dates, file formats, and grant numbers are used to ensure compliance with both internal and external data-sharing policies.

## **Why You Should Contribute DataDSP Entries**
Contributing DataDSP entries benefits your research and projects by improving data accessibility, organization, and compliance. With complete entries, you enhance collaboration opportunities, simplify data sharing, and increase the impact and visibility of your datasets in research communities. By ensuring your data is properly documented and discoverable, you also reduce administrative burdens during audits, grant reporting, and data requests.

### **Who Should Be Contributing DataDSP Entries?**
1. **Principal Investigators (PIs)** – Gain recognition for your research by making your data easily accessible and well-documented, improving citation potential and collaboration opportunities.  
2. **Data Managers** – Ensure efficient data organization and compliance, minimizing time spent addressing data queries or audits.  
3. **Research Coordinators** – Streamline project workflows by contributing accurate metadata, reducing delays in data sharing and project reporting.  
4. **Consortium Members** – Enhance collaboration by contributing standardized data entries, ensuring that datasets are usable across multiple institutions and research projects. 


## Download Template
Download the [DataDSP CSV template](https://github.com/mc2-center/data-models/raw/main/templates/DataDSP.csv) for streamlined data entry, ensuring that all required fields are filled out.


## Example Data Entry
The table below includes sample values to demonstrate proper attribute usage.

| **Attribute**               | **Example Value**                                                                                                      |
|-----------------------------|-------------------------------------------------------------------------------------------------------------------------|
| DSP Dataset Name            | DSP_Dataset_Lung_Research_2021                                                                                          |
| DSP Dataset Alias           | Syn123456                                                                                                          |
| DSP Dataset Assay           | 3D Bioprinting                                                                                                           |
| DSP Dataset Species         | Asian Elephant                                                                                                           |
| DSP Dataset File Formats    | CSV, JSON                                                                                                                |
| DSP Planned Upload Date     | 2022-12-01                                                                                                               |
| DSP Dataset Grant Number    | CA209971                                                                                                                  |
| DSP Dataset Description     | "A quick description of your data sharing plan."            |
| DSP Dataset Destination     | "/home/user/datasets/dsp_output"                                                                                          |
| DSP Dataset Url             | https://www.example.com/dataset/dsp1234                                                                                   |


## Full Field Reference

Below is the full field reference table with attributes and their descriptions.


[⤓ Download template](https://github.com/mc2-center/data-models/raw/main/templates/DataDSP.csv)

{{ read_csv('sharingPlans/template.csv') }}