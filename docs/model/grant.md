A **Grant** entry captures essential details about funded research projects. Grants support various scientific studies, often involving multiple institutions, investigators, and funding agencies. Accurate grant records facilitate resource tracking, collaboration, and compliance with funding requirements.

This section outlines how to create a valid grant entry by defining key metadata, including grant type, institution, investigator, and project duration. Proper grant documentation ensures projects can be easily identified, referenced, and managed across research initiatives.

### **Why You Should Contribute Grant Entries**
Contributing grant entries ensures that funding sources, research objectives, and project details are transparent and traceable. This helps stakeholders monitor the progress of research initiatives, improves collaboration by connecting teams with similar funding, and supports compliance with reporting and accountability requirements. Accurate grant documentation also strengthens your ability to secure future funding by showcasing project outcomes and impacts.

#### **Who Should Be Contributing Grant Entries?**
1. **Principal Investigators (PIs)** – Document your funded projects to highlight ongoing research, objectives, and collaborations for visibility and accountability.  
2. **Research Administrators** – Maintain accurate grant records for compliance, progress monitoring, and reporting purposes.  
3. **Data Managers** – Link datasets, publications, and other project outputs to their corresponding funding sources to improve data traceability.  
4. **Grant Writers and Project Coordinators** – Provide detailed grant descriptions and metadata to showcase successful projects and support future funding proposals.  
5. **Consortia and Collaborative Project Leads** – Share grant details to foster transparency and facilitate multi-institutional research efforts across teams and disciplines.  


## Download Template
Use the [grant entry template](https://github.com/mc2-center/data-models/raw/main/templates/GrantView.csv) to streamline your data entry process. The template contains pre-defined required fields.

## Example Data Entry
The table below includes sample values to demonstrate proper attribute usage.

| **Attribute**              | **Example Value**                                                                                           |
|----------------------------|-------------------------------------------------------------------------------------------------------------|
| Grant Name                 | Research Grant for Cancer Genomics                                                                           |
| Grant Number               | CA209971                                                                                                     |
| Grant Abstract             | This grant funds research into genetic mutations contributing to aggressive cancers, focusing on tumor sequencing and gene expression analysis. |
| Grant Type                 | R37                                                                                                          |
| Grant View                 | Public                                                                                                       |
| Grant Theme Name           | Cancer Genomics, Precision Medicine                                                                          |
| Grant Institution Name      | Stanford University                                                                                         |
| Grant Institution Alias    | Stanford                                                                                                     |
| Grant Investigator         | Dr. John Smith                                                                                               |
| Grant Consortium Name      | HTAN                                                                                                         |
| GrantView_id               | GV12345                                                                                                      |
| Grant Synapse Team         | Team: Project Management, Permission: Edit                                                                   |
| Grant Synapse Project      | Synapse_ID: syn123456789, Grant_Name: NIH Brain Initiative Grant                                             |
| Grant Start Date           | 2022-01-01                                                                                                   |
| NIH RePORTER Link          | [NIH Reporter Project Info](https://projectreporter.nih.gov/project_info_description.cfm?aid=9813521)         |
| Duration of Funding        | 5 years                                                                                                      |
| Embargo End Date           | 2023-12-31                                                                                                   |



## Full Field Reference

Below is the full field reference table with attributes and their descriptions.

[⤓ Download template](https://github.com/mc2-center/data-models/raw/main/templates/GrantView.csv)

{{ read_csv('grant/template.csv') }}
