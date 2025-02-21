A **Publication** entry documents key information about scientific articles, research studies, and related materials. Properly recording publications ensures that data is easily traceable, searchable, and meets compliance standards. This section explains how to create a valid publication entry using required fields, templates, and example data for guidance.

### **Why You Should Contribute Publication Entries**
Contributing publication entries helps highlight research achievements, supports compliance with funding and reporting requirements, and ensures the discoverability of research outputs. Accurate publication records promote transparency, facilitate citation tracking, and showcase the impact of research to stakeholders. Including comprehensive details such as PubMed IDs, keywords, and grant numbers enables efficient data organization and cross-referencing within research ecosystems.

#### **Who Should Be Contributing Publication Entries?**
1. **Researchers and Authors** – Showcase your research contributions, improve visibility, and enable easier access to your published work through linked identifiers and citations.  
2. **Project Leads and Investigators** – Ensure that key publications resulting from funded research projects are documented, supporting compliance and progress reporting.  
3. **Research Administrators** – Maintain an accurate record of research outputs tied to grants and projects, which aids in tracking productivity and measuring impact.  
4. **Consortium Members** – Document publications across collaborations to highlight the collective research contributions and foster recognition within the scientific community.  
5. **Data Managers** – Link datasets, tools, and studies to their associated publications, providing full context for data provenance and enhancing traceability.  


## Download Template
Download the [publication entry template](https://github.com/mc2-center/data-models/raw/main/templates/PublicationView.csv) to streamline data entry and ensure required fields are completed.

## Example Data Entry
The table below includes sample values to demonstrate proper attribute usage.

| **Attribute**              | **Example Value**                                                                                                  |
|----------------------------|---------------------------------------------------------------------------------------------------------------------|
| Publication Doi            | [10.1016/j.cell.2016.02.013](https://doi.org/10.1016/j.cell.2016.02.013)                                             |
| Publication Journal        | *Nature Medicine*                                                                                                   |
| Pubmed Id                  | 12345678                                                                                                             |
| Pubmed Url                 | [PubMed Link](https://www.ncbi.nlm.nih.gov/pubmed/12345678)                                                          |
| Publication Title          | Role of Inflammatory Pathways in Cancer Progression                                                                 |
| Publication Year           | 2023                                                                                                                 |
| Publication Keywords       | Cancer Pathways, Inflammation, Tumor Microenvironment                                                               |
| Publication Authors        | Dr. Sarah Johnson, Dr. Michael Lee, Prof. Amanda Carter                                                              |
| Publication Abstract       | This study investigates inflammatory pathways and their impact on tumor progression and metastasis.                   |
| Publication Assay          | RNA Sequencing, In Vivo Bioluminescence                                                                              |
| Publication Tumor Type     | Lung Cancer, Breast Cancer                                                                                           |
| Publication Tissue         | Lung, Breast                                                                                                         |
| Publication Accessibility  | Open Access                                                                                                          |
| Publication View           | Online                                                                                                               |
| Publication Grant Number   | CA302123, CA301987                                                                                                    |
| Publication Dataset Alias  | GSE45678, DOI:10.1000/exampledataset                                                                                 |
| PublicationView_id         | PublicationView_45678                                                                                                |



## Full Field Reference

Below is the full field reference table with attributes and their descriptions.

[⤓ Download template](https://github.com/mc2-center/data-models/raw/main/templates/PublicationView.csv)

{{ read_csv('publication/template.csv') }}