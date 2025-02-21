A **Person** entry provides detailed information about individuals involved in research projects, including their affiliations, roles, and contributions. Accurate person records are essential for tracking collaborations, grant participation, working group involvement, and expertise within scientific initiatives.

This section outlines how to create a valid person entry, ensuring that all necessary details are captured to support identification, categorization, and data integration within  MC<sup>2</sup> Center-supported projects.

### **Why You Should Contribute Person Entries**
Contributing person entries enhances collaboration by providing clear visibility into expertise, roles, and project involvement within research teams. It helps connect researchers and contributors across projects, facilitates networking opportunities, and allows stakeholders to easily identify key personnel for collaboration or consultation. Maintaining up-to-date profiles also supports transparency and improves recognition of contributions to research initiatives.

#### **Who Should Be Contributing Person Entries?**
1. **Researchers and Scientists** – Document your roles, publications, and collaborations to increase visibility within your field and enable others to reach out for partnerships.  
2. **Project Leads and Investigators** – Ensure key contributors to your projects are credited and easily searchable for collaboration and reporting purposes.  
3. **Research Administrators** – Keep records of team members, roles, and working group participation to manage projects efficiently and facilitate compliance reporting.  
4. **Consortium Members** – Share your expertise and involvement across multiple projects or institutions to encourage cross-disciplinary collaboration and integration of research efforts.  
5. **Data Managers** – Link researchers to datasets, tools, and resources they have contributed to ensure traceability and acknowledgment of contributions.  

## Download Template
Download the [person entry template](https://github.com/mc2-center/data-models/raw/main/templates/PersonView.csv) for streamlined data entry, ensuring that all required fields are filled out.

## Example Data Entry
The table below includes sample values to demonstrate proper attribute usage.

| **Attribute**               | **Example Value**                                                                                     |
|-----------------------------|--------------------------------------------------------------------------------------------------------|
| Name                        | Dr. Jane Smith                                                                                        |
| Alternative Names           | Jane Doe, J. Smith                                                                                    |
| Email                       | janesmith@example.com                                                                                 |
| Url                         | [Professional Profile](https://www.example.com/janesmith)                                             |
| Orcid Id                    | 0000-0002-1825-0097                                                                                   |
| Synapse Profile Id          | SP98765432                                                                                            |
| Last Known Institution      | Stanford University                                                                                   |
| Working Group Participation | Cancer Metabolism, Genomic Data Integration                                                           |
| Chair Roles                 | Steering Committee, Ethics Review Panel                                                               |
| Consent For Portal Display  | Yes                                                                                                    |
| Portal Display              | TRUE                                                                                                   |
| Person View                 | Public Profile View                                                                                    |
| Person Grant Number         | CA202177, CA304599                                                                                     |
| Person Consortium Name      | HTAN, PDMC                                                                                            |
| Person Publications         | 25700473, 31245678                                                                                     |
| Person Datasets             | GSE12345, DOI:10.1000/sampledataset                                                                   |
| Person Tools                | R Studio, Python, Jupyter Notebook                                                                    |
| Person Educational Resources | Advanced Cancer Research Training, NIH Data Science Workshop                                           |
| PersonView_id               | PersonView_12345                                                                                       |

## Full Field Reference

Below is the full field reference table with attributes and their descriptions.

[⤓ Download template](https://github.com/mc2-center/data-models/raw/main/templates/PersonView.csv)

{{ read_csv('person/template.csv') }}