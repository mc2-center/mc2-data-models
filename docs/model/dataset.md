A **Dataset** refers to a structured collection of data that is organized for analysis, sharing, and research. Datasets can include various types of information such as demographic statistics, experimental results, or survey responses. In  MC<sup>2</sup> Center-supported projects, dataset entries help ensure data is properly categorized, easily retrievable, and compliant with sharing and storage requirements.

This model outlines key attributes that describe and manage datasets, including metadata about the data type, format, number of samples, and related grant information. By maintaining these attributes, datasets can be efficiently tracked and referenced within data repositories.

## **Why You Should Contribute Dataset Entries**
Contributing dataset entries helps ensure your data is easily discoverable, accurately documented, and compliant with research standards. By providing well-structured dataset metadata, you enhance opportunities for collaboration, boost data citation potential, and reduce administrative overhead for reporting and compliance. Well-documented datasets also help other researchers and stakeholders use your data effectively, increasing its long-term impact.

### **Who Should Be Contributing Dataset Entries?**
1. **Principal Investigators (PIs)** – Increase the visibility and impact of your research by contributing properly cataloged datasets, making it easier for others to cite and use your work.  
2. **Data Managers** – Improve data organization and retrieval, reducing time spent on requests for data clarification and documentation during audits.  
3. **Research Staff** – Simplify project reporting by ensuring that datasets are complete with accurate descriptions, grant associations, and metadata.  
4. **Collaborators and Partners** – Enhance data interoperability across multiple institutions by contributing standardized entries that support shared research initiatives. 


## Download Template
You can download the [dataset entry template](https://github.com/mc2-center/data-models/raw/main/templates/DatasetView.csv), which includes all required fields, to streamline the data entry process.

## Example Data Entry
The table below includes sample values to demonstrate proper attribute usage.

| **Attribute**           | **Example Value**                                                                                       |
|-------------------------|---------------------------------------------------------------------------------------------------------|
| Dataset Name            | RNA Sequencing of Lung Cancer Samples 2021                                                                                   |
| Dataset Alias           | GSE56789                                                                                             |
| Dataset Description     | This dataset contains RNA sequencing data from 200 lung cancer samples, including gene expression profiles and patient clinical data. It is designed to study differential gene expression and mutation burden across tumor stages.             |
| Dataset Url             | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56789                                    |
| Dataset Assay           | RNA Sequencing                                                                                                    |
| Dataset Species         | Homo sapiens                                                                                                   |
| Dataset Tumor Type      | Glioblastoma                                                                                          |
| Dataset Tissue          | Lung                                                                                          |
| Dataset File Formats    | CSV, PDF                                                                                                |
| Dataset Grant Number    | CA209971                                                                                                |
| Dataset Pubmed Id       | Not applicable                                                                                          |
| Dataset View            | Table                                                                                                   |
| DatasetView_id          | DatasetView_12345                                                                                       |



## Full Field Reference

Below is the full field reference table with attributes and their descriptions.

[⤓ Download template](https://github.com/mc2-center/data-models/raw/main/templates/DatasetView.csv)

{{ read_csv('dataset/template.csv') }}
