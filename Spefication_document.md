

## 

## 

# **Specification Document**

# 

# InfectioGIT : A Structured Repository of Models and Metadata for Human Infectious Diseases

**Students : Aya Sebbah \- Diana Bravais \- Thalie Holmiere**  
**Master 1 – Bioinformatics & Systems Biology**   
**Supervisors : Pr. Anna Niarakis & Dr Virginie Jouffret**  
**Affiliation : CBI**

**2025-2026**

# Table of Content

[**Table of Content	1**](#table-of-content)

[**I. Project Context	2**](#project-context)

[**II. General Objective	3**](#general-objective)

[**III. Key Steps	4**](#key-steps)

[1\. Development of a curated catalogue of infectious diseases	4](#development-of-a-curated-catalogue-of-infectious-diseases)

[A. Selection Criteria	4](#selection-criteria)

[a. Public health priority selection	4](#public-health-priority-selection)

[b. Algorithmic prioritization	4](#algorithmic-prioritization)

[B. Selected Diseases (initial panel)	4](#selected-diseases-\(initial-panel\))

[2\. Model Retrieval, Filtering, and Metadata Standardization Pipeline	5](#model-retrieval,-filtering,-and-metadata-standardization-pipeline)

[A. Definition of Search Keywords	6](#definition-of-search-keywords)

[B. Retrieval of Models from BioModels	6](#retrieval-of-models-from-biomodels)

[C. State of the Art: Disease-Specific Model Landscape	7](#state-of-the-art:-disease-specific-model-landscape)

[a. Quantitative Assessment (Model Volume)	7](#quantitative-assessment-\(model-volume\))

[b. Model Filtering and Selection	7](#model-filtering-and-selection)

[c. Identification of Knowledge Gaps	7](#identification-of-knowledge-gaps)

[D. Metadata Acquisition Structuring and Standardization	8](#metadata-acquisition-structuring-and-standardization)

[A unified metadata schema will be developed based on :	8](#a-unified-metadata-schema-will-be-developed-based-on-:)

[E. Automated Metadata Enrichment	8](#automated-metadata-enrichment)

[**IV. Technical Materials and Architecture	9**](#technical-materials-and-architecture)

[1\. Data Architecture and Repository Structure	9](#data-architecture-and-repository-structure)

[A. Materials	9](#materials)

[B. InfectioGIT GitHub Repository Organization	9](#infectiogit-github-repository-organization)

[C. Technical Components and Material	10](#technical-components-and-material)

[**V. Constraints and Risks	10**](#constraints-and-risks)

[1\. Functional Constraints	10](#functional-constraints)

[2\. Organizational Constraints	10](#organizational-constraints)

[3\. Technical and Data Constraints	11](#technical-and-data-constraints)

[**VI. Final Deliverables	11**](#final-deliverables)

[**Bibliography	12**](#bibliography)

[**Annexes	14**](#annexes)

## 

1. # Project Context

   

InfectioGIT lies at the interface between computational biology and infectious disease research. It is developed as a sub-project within a broader Digital Twins initiative \[1\], whose long-term ambition is to build dynamic virtual representations of host–pathogen interactions and immune responses. By combining mechanistic knowledge with computational modelling, infectious-disease digital twins could help us understand, simulate, and ultimately predict how infections emerge, evolve within the host, and respond to interventions.

This project was selected in the context of EUR UNITEID, a program built on the One Health \[2\] principle, which considers human health, animal health, and environmental factors as tightly interconnected. Infectious diseases are a major One Health challenge because many pathogens circulate through complex ecosystems involving reservoirs, vectors, and environmental drivers. In particular, West Nile virus \[3\], which is the shared thematic focus of the UNITEID transdisciplinary initiative, illustrates this perfectly by involving  bird reservoirs, mosquito vectors, incidental hosts such as humans and horses, and transmission dynamics strongly shaped by seasonality, climate, and habitat. Addressing such systems requires not only biological expertise, but also contributions from data science, modelling, ecology, public health, and digital infrastructure design making the effort inherently transdisciplinary across multiple Master’s programs. 

A key obstacle to building infectious-disease digital twins is that relevant computational models remain scattered across publications and repositories, described with inconsistent metadata[^1], and difficult to compare, reproduce, or reuse. Yet these models are essential: they provide mechanistic formalisms to represent processes across scales (molecular, cellular, tissue, host, and sometimes population levels), and can be encoded in standard formats such as SBML \[4\] for quantitative biochemical networks or SBML-qual for qualitative and logical representations of regulatory and signalling systems. Without a structured way to organise these resources, the field loses time reinventing models rather than integrating and extending them.

Beyond long-term research goals, a centralised and well-annotated repository becomes especially critical during emerging outbreaks and health emergencies. In crisis situations, time is a limiting factor : researchers and decision-makers need rapid access to trustworthy models and structured knowledge to understand how a pathogen affects the host, anticipate disease dynamics, and evaluate potential interventions. By providing a single resource where models, metadata, and documentation are gathered in a consistent way. It overcomes the limited reusability of dispersed data, it saves researchers time in model discovery, comparison, and reuse.

To address this gap, InfectioGIT is designed as an open, FAIR\[5\]\-compliant, and community-driven resource that gathers, annotates, and organises computational models related to infectious diseases and host–pathogen interactions. It merges a GitHub[^2] repository with a lightweight structured database, creating a unified ecosystem for model discovery, standardised annotation, interoperability, and reuse. By centralising these modelling building blocks and making them easier to curate collaboratively, InfectioGIT aims to provide the scalable digital foundation needed to support infectious-disease digital twins, including One Health–oriented applications such as West Nile virus modelling within the EUR UNITEID framework.

**Host team :** 

This project is supervised by Professor Anna Niarakis and Dr. Virginie Jouffret at the Center for Integrative Biology (CBI) in Toulouse. Professor Niarakis is a specialist in systems biology and leads the Computational Systems Biology team. Dr. Jouffret is an expert in bioinformatics and biostatistics from the BigA platform. Their collaboration at CBI provides the necessary expertise in both biological modeling and data processing for this project.

The CBI is an international multidisciplinary research institute dedicated to understanding biological systems through integrative and multi-scale approaches. It federates several laboratories, including the host laboratory Molecular, Cellular and Developmental Biology (MCD), and addresses major health-related societal challenges by combining experimental biology with quantitative and computational methods.

Within this scientific environment, the host team contributes to the development of integrative computational approaches at the crossroads of mathematics, computer science, and bioinformatics. Their expertise covers the design of reproducible workflows, the structuring of biological knowledge, and the development of tools that make computational resources easier to access, assess, and reuse by the community.

2. # General Objective

The objective of the InfectioGIT project is to collect, organize, standardize, and make accessible computational models describing human infectious diseases, using biomedical ontologies and the BioModels[^3] \[6\] repository as the primary source of mechanistic models.

The project aims to build a centralized, structured, and queryable resource of SBML models relevant to emerging infectious diseases.

3. # Key Steps

1. ## **Development of a curated catalogue of infectious diseases**

1. ### **Selection Criteria**

Diseases will be selected based on two complementary approaches:

1. #### **Public health priority selection**

Using **emergent diseases** which refers to diseases whose incidence or geographic range is increasing, that can (re-)appear unexpectedly, and/or that pose a growing outbreak risk due to factors like climate-driven vector expansion, animal–human spillover, travel, and gaps in population immunity, relevant for arboviruses (dengue, chikungunya, West Nile) and other outbreak-prone infections.

- WHO[^4] global health priorities  
- National priorities[^5] defined by Santé Publique France[^6]  
- Pandemic preparedness frameworks

2. #### **Algorithmic prioritization** 

Through **high-modeling diseases**, those with an unusually large number of published computational models and rich datasets, often because they have had major outbreaks, sustained surveillance, and strong research investment (e.g., COVID-19, dengue). Including them is useful for benchmarking pipelines and comparing model types at scale.

##### 

A Python-based scoring algorithm is to be developed to rank diseases based on:

- Environmental exposure factors  
- Host diversity  
- Scarcity of existing computational models  
- Urgency of epidemiological reporting

2. ### **Selected Diseases (initial panel)** 

National (France) : 

- **West Nile virus** : a mosquito-borne flavivirus that can cause fever and, in some cases, severe neurological disease such as meningitis or encephalitis. Most infections are asymptomatic.  
- **Dengue** \[7\] : a mosquito-borne viral disease (dengue virus) causing high fever, severe headache, and muscle/joint pain; some cases progress to severe dengue with bleeding and shock.  
- **Chikungunya** \[8\] : a mosquito-borne alphavirus infection characterized by sudden fever and intense joint pain, which can persist for weeks to months. It may also cause rash and fatigue.

International: 

- **Mpox** \[9\] : a zoonotic orthopoxvirus infection that causes fever and a distinctive rash with lesions, often with swollen lymph nodes. Transmission occurs through close contact with infected people or animals and contaminated materials.  
- **Avian influenza** (bird flu) and **Measles Influenza** \[10\] : a virus primarily circulating in birds that can sometimes infect humans, ranging from mild illness to severe respiratory disease. Human cases are typically linked to exposure to infected birds or contaminated environments.  
- **COVID-19** \[11\] (optional inclusion for data volume and benchmarking) : respiratory infectious disease caused by SARS-CoV-2, with symptoms ranging from mild upper-respiratory illness to severe pneumonia and multi-organ complications. It spreads mainly via respiratory droplets/aerosols.

2. ## **Model Retrieval, Filtering, and Metadata Standardization Pipeline**

To extract data, we will primarily prioritize resources from the **Disease Ontology \[12\]** (DO), a standardized hierarchy of human diseases that provides unique identifiers to eliminate naming ambiguity and enable precise disease mapping. It serves as the project's semantic backbone for the consistent classification and indexing of computational models. For each disease, we will collect two distinct ontology identifiers: **DOID[^7]** and **MONDO[^8]** \[13\]. The Disease Ontology identifiers (**DOID**) will be sourced from the HumanDO.obo file, which a custom Python script parses to extract identifiers corresponding to each disease of interest. Simultaneously, we will integrate **MONDO** identifiers, a unified disease ontology that harmonizes multiple resources into a single, cohesive hierarchy to improve cross-species data integration. This dual-ontology strategy is essential for ensuring maximal model coverage within the BioModels database, as entries are inconsistently annotated, some relying exclusively on DOID while others use MONDO. By retaining both, we will avoid missing relevant entries and ensure a comprehensive and high-quality final database.

   /IDs/  
   dengue\_DOID.txt  
   dengue\_MONDO.txt  
   chikungunya\_DOID.txt  
   chikungunya\_MONDO.txt  
Figure 1 : Example of the file structure for extracted data

1. ### **Definition of Search Keywords**  

**For each disease, a list of standardized search terms will be defined, including :**

- Synonyms  
- Scientific names  
- Alternative spellings  
- Acronyms

This step is essential because BioModels annotations may vary in naming conventions.

2. ### **Retrieval of Models from BioModels** 

**Data access** : Models will be retrieved using :

- BioModels REST API  
- Python Bioservices[^9] package

**Retrieval Pipeline** : A Python script will be developed to perform  the following:

1. Queries BioModels using:  
- DOID identifiers  
- MONDO identifiers  
- keyword list  
2. Filter relevant models (based on number of citations of the linked article)  
3. Download SBML files  
   /models/dengue/  
      BIOMDxxxx.xml  
      metadata\_dengue.json  
      model\_list.txt

Figure 2 : Example of output files for each disease

3. ### **State of the Art: Disease-Specific Model Landscape** {#state-of-the-art:-disease-specific-model-landscape}

A systematic state-of-the-art review will be performed to refine the selection of **30 to 60 models**. This phase will evaluate the modeling landscape based on three pillars :

1. #### **Quantitative Assessment (Model Volume)** 

**Quantification of the available SBML/ODE models** in BioModels for each target disease. It will be done to identify High-Modeling Diseases (e.g., COVID-19, Dengue) to serve as pipeline benchmarks.

**Redundancy Management** : If an excess of models is found, a selection will be made based on scientific relevance (citation number on articles) and to select a diverse array of computational models approaches.

2. #### **Model Filtering and Selection** 

The project will evaluate the **scientific weight** of candidate models by cross-referencing BioModels entries with **PubMed/Google** Scholar **citation counts**.

**Selection Criteria** : Models are filtered based on:

- Biological relevance  
  - Citation-Based Filtering : A Python script using Biopython[^10]’s Entrez retrieves publication data and filters models based on:  
    - Number of citations  
    - Journal quality  
    - Publication year  
- Annotation completeness  
- Availability of an associated scientific publication

Priority will be given to "standard" models (high citation rates), ensuring the catalogue is built on validated and widely recognized mathematical frameworks.

3. #### **Identification of Knowledge Gaps** 

A key goal will be to identify and document **Modeling Gaps** The project will specifically look for:

- **Scarcity Gaps** : Diseases with high outbreak risk but few computational models (e.g., Mpox, West Nile).  
- **Biological Gaps** : Lack of models integrating specific factors like host diversity or environmental exposure.

4. ### **Metadata Acquisition Structuring and Standardization** 

**Metadata Sourcing** :

- BioModels  metadata  
- COMBINE archive[^11] [\[14\]](https://sciwheel.com/work/citation?ids=18676809&pre=&suf=&sa=0) metadata   
- Disease Ontology (DOID)  
- MONDO ontology

**Metadata Standard** : 
### A unified metadata schema will be developed based on :

- COMBINE Archive metadata standards  
- JSON[^12]  structured format  
- FAIR[^13] data principles  
  {  
    "model\_id": "BIOMD0000000957",  
    "disease": "Dengue",  
    "ontology\_ids": \["DOID:12365", "MONDO:0001234"\],  
    "model\_type": "SIR",  
    "organism": "Human",  
    "publication": "...",  
    "authors": "...",  
    "year": 2020  
  }

Figure 3 : Example of Metadata Structure

5. ### **Automated Metadata Enrichment**

To address the lack of standardized annotations in non-curated models, the project will implement an automated enrichment pipeline.

**Enrichment Process** :

1. **Gap Identification** : A Python script will detect models with missing or incomplete ontology links (DOID/MONDO).

2. **Automated Injection** : The script will parse source files from the Disease Ontology and MONDO to retrieve full hierarchical data (synonyms, parent terms, and cross-references).

3. **Metadata Update** : These attributes will be automatically injected into the local metadata.json files of each model.

This will ensure that all models, regardless of their original curation status, reach a uniform standard of quality and remain fully searchable within the InfectioGIT ecosystem.

4. # Technical Materials and Architecture

1. ## **Data Architecture and Repository Structure** 

The primary goal will be to establish a centralized, open-source **InfectioGIT GitHub repository**. This infrastructure will be designed to host a standardized collection of infectious disease models, ensuring they are logically organized, machine-readable, and fully compliant with **FAIR data principles**. By partitioning models, metadata, and ontologies, the repository will provide a scalable environment for automated biological research.

1. ### **Technical Stack & Components** 

- **Python 3.13**, **Bioservices** and **BioPython:** Core language and specialized bioinformatics library.  
- **SQLite:** Lightweight SQL database.  
- **AI Assistants (GPT-4o & Gemini 2.0 Flash/Pro):** Used for code optimization and technical drafting.

**Development Environment and Language Rationale**

**Programming Language:** The entire data pipeline will be developed in **Python 3.14.3**.

Python was selected as the core language due to the team’s advanced proficiency and its offers of an extensive ecosystem of specialized **bioinformatics packages** (e.g., *Biopython*, *Bioservices*) that are essential for interfacing with biological databases and processing XML/SBML formats.

- **SBML Models** (.xml) : The scientific core of the project. These machine-readable files provide the mechanistic and mathematical descriptions (ODEs, reaction networks) of infectious processes, ensuring simulation compatibility and reproducibility.

- **Metadata Files** (.json) : Descriptive files that facilitate efficient indexing, searching, and filtering. They are essential for FAIR compliance, allowing the repository to interact with external bioinformatics services without parsing complex XML code.

- **Ontology Files** (.obo / .owl) : Semantic resources from DOID and MONDO. They provide a standardized vocabulary to resolve naming ambiguities, support hierarchical classification, and enable precise disease mapping across global databases.

- **Python Scripts** (.py) : Automation tools used to manage the data lifecycle. These scripts handle model retrieval from BioModels, metadata enrichment, and structural quality control, ensuring the scalability of the curation workflow.

- **SQLite Database** : A lightweight, portable RDBMS used to index metadata. It was selected for its low deployment overhead, providing a fast and queryable index of the entire model collection.

2. ### **InfectioGIT GitHub Repository Organization** 

The repository is structured as follows :

/InfectioGIT  
   /disease\_name/ (e.g., /dengue/): contains disease-specific data.  
       /models/  
	SBML files (mechanistic and mathematical descriptions)  
       /metadata/  
JSON files (indexing, searchability, and interoperability)  
   /ontologies/ for standardized disease mapping.  
       DOID.obo  
       MONDO.owl  
   /scripts/ Python automation tools for retrieval and enrichment.  
   /docs/ Technical documentation and user guides.  
Figure 4 : Organizational structure of the InfectioGIT repository

4. # Constraints and Risks 

1. ##  **Functional Constraints**  

The technical framework of the project is strictly defined by the use of Python for automation and data processing scripts, alongside SQLite for database management. These tools are selected to guarantee interoperability and efficiency when handling SBML and JSON file formats, ensuring the resource remains accessible and easy to maintain.

2. ##  **Organizational Constraints**  

Team coordination is primarily shaped by our dual enrollment in distinct academic tracks and the simultaneous management of the EUR UNITEID program alongside a second tutored project. While the five-month timeframe for this tutored project is ambitious for a project of this scope, we have leveraged these constraints to foster a highly efficient and disciplined workflow. We recognize that our commitment to two academic paths naturally limits the time available compared to a standard curriculum, which has in turn reinforced our need for impeccable organization. Consequently, we rely on collaborative platforms like GitHub to maximize our collective output and ensure steady progress despite the condensed schedule and heavy academic workload.

3. ##  **Technical and Data Constraints**

The availability of information is a critical constraint, as BioModels does not cover every infectious virus of interest for this study. This necessitates the documentation and integration of alternative resources such as NCBI, a comprehensive public database for genomic and biomedical information, Zenodo, a multi-disciplinary open-access repository for research datasets, or PubMed, a vast scholarly search engine for academic literature.

Currently, BioModel is open-source but may no longer be available because it is based in the United States due to geopolitical issues.

6. # Final Deliverables 

- A curated collection of SBML infectious disease models and metadata  
- A public GitHub repository (InfectioGIT)  
- A structured and queryable database with a full technical documentation to help researchers to access data   
- A specifications document, this document outlines the technical specifications, functional requirements, and overall framework of the project.  
- A reproducible pipeline for other diseases (i.e., a standardized workflow that can be reused and adapted for additional diseases).

# Bibliography 

1\. Niarakis, A., Laubenbacher, R., An, G., Ilan, Y., Fisher, J., Flobak, Å., Reiche, K., Rodríguez Martínez, M., Geris, L., Ladeira, L., Veschini, L., Blinov, M. L., Messina, F., Fonseca, L. L., Ferreira, S., Montagud, A., Noël, V., Marku, M., Tsirvouli, E., … Glazier, J. A. (2024). Immune digital twins for complex human pathologies: applications, limitations, and challenges. NPJ Systems Biology and Applications, 10(1), 141\. https://doi.org/10.1038/s41540-024-00450-5  
2\. Sinclair, J. R. (2019). Importance of a One Health approach in advancing global health security and the Sustainable Development Goals. Revue Scientifique et Technique (International Office of Epizootics), 38(1), 145–154. https://doi.org/10.20506/rst.38.1.2949  
3\. Martin, M.-F., & Simonin, Y. (2019). \[West Nile virus historical progression in Europe\]. Virologie (Montrouge, France), 23(5), 265–270. https://doi.org/10.1684/vir.2019.0787  
4\. Chaouiya, C., Bérenguier, D., Keating, S. M., Naldi, A., van Iersel, M. P., Rodriguez, N., Dräger, A., Büchel, F., Cokelaer, T., Kowal, B., Wicks, B., Gonçalves, E., Dorier, J., Page, M., Monteiro, P. T., von Kamp, A., Xenarios, I., de Jong, H., Hucka, M., … Helikar, T. (2013). SBML qualitative models: a model representation format and infrastructure to foster interactions between qualitative modelling formalisms and tools. BMC Systems Biology, 7, 135\. https://doi.org/10.1186/1752-0509-7-135  
5\. Wilkinson, M. D., Dumontier, M., Aalbersberg, I. J. J., Appleton, G., Axton, M., Baak, A., Blomberg, N., Boiten, J.-W., da Silva Santos, L. B., Bourne, P. E., Bouwman, J., Brookes, A. J., Clark, T., Crosas, M., Dillo, I., Dumon, O., Edmunds, S., Evelo, C. T., Finkers, R., … Mons, B. (2016). The FAIR Guiding Principles for scientific data management and stewardship. Scientific Data, 3, 160018\. https://doi.org/10.1038/sdata.2016.18  
6\. Le Novère, N., Bornstein, B., Broicher, A., Courtot, M., Donizelli, M., Dharuri, H., Li, L., Sauro, H., Schilstra, M., Shapiro, B., Snoep, J. L., & Hucka, M. (2006). BioModels Database: a free, centralized database of curated, published, quantitative kinetic models of biochemical and cellular systems. Nucleic Acids Research, 34(Database issue), D689-91. https://doi.org/10.1093/nar/gkj092  
7\. Guzman, M. G., & Harris, E. (2015). Dengue. The Lancet, 385(9966), 453–465. https://doi.org/10.1016/S0140-6736(14)60572-9  
8\. Maure, C., Khazhidinov, K., Kang, H., Auzenbergs, M., Moyersoen, P., Abbas, K., Santos, G. M. L., Medina, L. M. H., Wartel, T. A., Kim, J. H., Clemens, J., & Sahastrabuddhe, S. (2024). Chikungunya vaccine development, challenges, and pathway toward public health impact. Vaccine, 42(26), 126483\. https://doi.org/10.1016/j.vaccine.2024.126483  
9\. Hou, W., Wu, N., Liu, Y., Tang, Y., Quan, Q., Luo, Y., & Jin, C. (2025). Mpox: Global epidemic situation and countermeasures. Virulence, 16(1), 2457958\. https://doi.org/10.1080/21505594.2025.2457958  
10\. Bi, Y., Yang, J., Wang, L., Ran, L., & Gao, G. F. (2024). Ecology and evolution of avian influenza viruses. Current Biology, 34(15), R716–R721. https://doi.org/10.1016/j.cub.2024.05.053  
11\. Ochani, R., Asad, A., Yasmin, F., Shaikh, S., Khalid, H., Batra, S., Sohail, M. R., Mahmood, S. F., Ochani, R., Hussham Arshad, M., Kumar, A., & Surani, S. (2021). COVID-19 pandemic: from origins to outcomes. A comprehensive review of viral pathogenesis, clinical manifestations, diagnostic evaluation, and management. Le Infezioni in Medicina : Rivista Periodica Di Eziologia, Epidemiologia, Diagnostica, Clinica e Terapia Delle Patologie Infettive, 29(1), 20–36.  
12\. Schriml, L. M., Mitraka, E., Munro, J., Tauber, B., Schor, M., Nickle, L., Felix, V., Jeng, L., Bearer, C., Lichenstein, R., Bisordi, K., Campion, N., Hyman, B., Kurland, D., Oates, C. P., Kibbey, S., Sreekumar, P., Le, C., Giglio, M., & Greene, C. (2019). Human Disease Ontology 2018 update: classification, content and workflow expansion. Nucleic Acids Research, 47(D1), D955–D962. https://doi.org/10.1093/nar/gky1032  
13\. Vasilevsky, N. A., Toro, S., Matentzoglu, N., Flack, J. E., Mullen, K. R., Hegde, H., Gehrke, S., Whetzel, P. L., Shwetar, Y., Harris, N. L., Ngu, M. S., Alyea, G. L., Kane, M. S., Roncaglia, P., Sid, E., Thaxton, C. L., Wood, V., Abraham, R. S., Achatz, M. I., … Haendel, M. A. (2025). Mondo: integrating disease terminology across communities. Genetics. https://doi.org/10.1093/genetics/iyaf215  
14\. Gennari, J. H., König, M., Misirli, G., Neal, M. L., Nickerson, D. P., & Waltemath, D. (2021). OMEX metadata specification (version 1.2). Journal of Integrative Bioinformatics, 18(3). https://doi.org/10.1515/jib-2021-0020

# Annexes 

Annexe 1 : Estimated Project Timeline and Key Deliverables

[^1]:  [Metadata :](https://www.dublincore.org/resources/metadata-basics/?utm_source=chatgpt.com) Metadata are structured information that describe a resource or dataset (e.g., title, author, date, format, provenance) to support discovery, understanding, and reuse. 

[^2]:  [GitHub:](https://docs.github.com/en/get-started/start-your-journey/about-github-and-git?utm_source=chatgpt.com)  a web-based platform for hosting **Git repositories** and supporting **collaborative software development** (version control, code review, issues, project management). 

[^3]:  [BioModels](https://www.ebi.ac.uk/training/online/courses/biomodels-quick-tour/what-is-biomodels/?utm_source=chatgpt.com) is a free, public repository that stores, shares, and provides access to published computational/mathematical models of biological and biomedical systems, with many models manually curated to ensure reproducibility and correct annotation.

[^4]:  [WHO](https://www.who.int/fr/news-room/fact-sheets/detail/the-top-10-causes-of-death) is the United Nations agency that connects nations, partners and people to promote health, keep the world safe and serve the vulnerable – so everyone, everywhere can attain the highest level of health. 

[^5]:  https://anrs.fr/recherche/projets-de-recherche/

[^6]:  [Santé Publique France](https://www.santepubliquefrance.fr/maladies-et-traumatismes/maladies-a-transmission-vectorielle/chikungunya/documents/bulletin-national/chikungunya-dengue-zika-et-west-nile-en-france-hexagonale.-bulletin-de-la-surveillance-renforcee-du-5-novembre-2025?utm_source=chatgpt.comhttps://www.santepubliquefrance.fr/maladies-et-traumatismes/maladies-a-transmission-vectorielle/chikungunya/documents/bulletin-national/chikungunya-dengue-zika-et-west-nile-en-france-hexagonale.-bulletin-de-la-surveillance-renforcee-du-5-novembre-2025?utm_source=chatgpt.com) is the national public health agency of France. It is responsible for monitoring, preventing, and controlling health risks across the country. The agency conducts research, collects health data, provides expertise to authorities, and promotes health education to improve population health. It plays a key role in responding to epidemics and emerging health threats.

[^7]: [**DOID (Disease Ontology):**](https://pubmed.ncbi.nlm.nih.gov/25348409/)  A stable identifier from the Disease Ontology (DO) used to uniquely label a human disease concept in a curated, hierarchical vocabulary (with definitions, synonyms, and cross-references).

[^8]: [**MONDO (Mondo Disease Ontology ID):**](https://mondo.readthedocs.io/en/latest/) A harmonized disease identifier from the Mondo Disease Ontology that integrates and maps disease concepts across multiple sources (e.g., DOID, Orphanet, OMIM), providing unified terms, synonyms, and cross-references for interoperability.

[^9]:  [BioServices](https://bioservices.readthedocs.io/en/main/?utm_source=chatgpt.com) is a Python framework that provides programmatic access to major bioinformatics web services (e.g., UniProt, KEGG, BioModels) and helps build wrappers for both REST and SOAP/WSDL services

[^10]:  [Biopython](https://biopython.org/) is a set of freely available tools for biological computation written in [Python](http://www.python.org/) by an international team of developers.

[^11]:  [COMBINE Archive](https://pmc.ncbi.nlm.nih.gov/articles/PMC4272562/?utm_source=chatgpt.com) is a single packaged file (OMEX, a ZIP container) that bundles all files needed to share and reproduce a computational modeling project (e.g., model files, simulation descriptions, data), together with a manifest and optional metadata.

[^12]:   JSON : JavaScript Object Notation is a lightweight, text-based data format for representing structured data using name–value pairs and ordered lists, designed to be easy for humans to read and write and for machines to parse and generate.

[^13]:  [FAIR data:](https://pmc.ncbi.nlm.nih.gov/articles/PMC4792175/?utm_source=chatgpt.com) Data (and associated metadata) that are managed so they are Findable, Accessible, Interoperable, and Reusable (FAIR), enabling reliable discovery, integration, and reuse by both humans and machines.