A **Study** entry documents key information about research projects, including methodologies, participants, and ethical considerations. Properly maintaining study records ensures that data is compliant with research standards, traceable, and usable across various scientific fields. This section explains how to create and maintain a valid study entry using the required fields, templates, and example data provided.


### **Why You Should Contribute Study Entries**
Contributing study entries ensures that critical research details such as methodology, investigators, and related datasets are documented and accessible. This fosters transparency, collaboration, and reproducibility in research. A well-maintained study entry also helps stakeholders understand the scope, purpose, and outcomes of research efforts, enabling better data integration and analysis across projects.

#### **Who Should Be Contributing Study Entries?**
1. **Principal Investigators (PIs)** – Provide clear documentation of your studies to enhance visibility and demonstrate research contributions for funding agencies and collaborators.

2. **Research Coordinators and Project Leads** – Ensure proper record-keeping of studies under your supervision, supporting compliance, reporting, and resource allocation.

3. **Data Managers** – Create study entries to link related datasets, publications, and tools, improving data provenance and facilitating cross-references within research ecosystems.

4. **Consortium Participants** – Contribute study information to showcase collaborative research efforts and track joint achievements.

5. **Funding Agencies and Grant Monitors** – Monitor study progress by ensuring that all funded research is accurately documented and accessible.

6. **Ethics and Compliance Officers** – Maintain study entries to verify that all necessary de-identification methods, IRB protocols, and data use policies are in place.


## Download Template
To streamline the process, download the [study entry template](https://github.com/mc2-center/data-models/raw/main/templates/Study.csv) for standardized data entry.


## Example Data Entry
The table below includes sample values to demonstrate proper attribute usage.

| **Attribute**                       | **Example Value**                                                                                                    |
|--------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| Study                                | Biology                                                                                                               |
| Study Name                           | Effects of Diet and Exercise on Obesity                                                                               |
| Study Description                    | Analysis of Cardiovascular Response during Exercise                                                                   |
| Study Investigator                   | Dr. Jane Doe, PhD in Nutrition Science, University of Texas                                                           |
| Study Reuse Statement                | Data from this study may be reused under conditions of proper citation and ethical approval.                           |
| Study_id                             | STUDY_2024_OBESITY_EXERCISE                                                                                           |
| Study Number of Participants         | 5000                                                                                                                  |
| Study De-identification Method Type  | Manual                                                                                                                |
| Study De-identification Method Description | Personal identifiers such as names and dates were removed, and randomization techniques were applied.               |
| Study De-identification Method Software | Safe Harbor Privacy Software                                                                                          |
| Study dbGaP Accession Id             | phs000424.v7.p2                                                                                                       |
| Study License                        | CC BY-NC 4.0                                                                                                          |
| Study Data Use Codes                 | IRB, PUB, HMB                                                                                                         |


## Full Field Reference

Below is the full field reference table with attributes and their descriptions.

[⤓ Download template](https://github.com/mc2-center/data-models/raw/main/templates/Study.csv)

{{ read_csv('study/template.csv') }}