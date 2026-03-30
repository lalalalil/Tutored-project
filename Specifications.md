# Specification Document

## InfectioGIT : A Structured Repository of Models and Metadata for Human Infectious Diseases

---

# Project Context

InfectioGIT lies at the interface between computational biology and infectious disease research. It is developed as a sub-project within a broader Digital Twins initiative [1], whose long-term ambition is to build dynamic virtual representations of host–pathogen interactions and immune responses. By combining mechanistic knowledge with computational modelling, infectious-disease digital twins could help us understand, simulate, and ultimately predict how infections emerge, evolve within the host, and respond to interventions.

This project was selected in the context of EUR UNITEID, a program built on the One Health [2] principle, which considers human health, animal health, and environmental factors as tightly interconnected. Infectious diseases are a major One Health challenge because many pathogens circulate through complex ecosystems involving reservoirs, vectors, and environmental drivers. In particular, West Nile virus [3], which is the shared thematic focus of the UNITEID transdisciplinary initiative, illustrates this perfectly by involving bird reservoirs, mosquito vectors, incidental hosts such as humans and horses, and transmission dynamics strongly shaped by seasonality, climate, and habitat. Addressing such systems requires not only biological expertise, but also contributions from data science, modelling, ecology, public health, and digital infrastructure design making the effort inherently transdisciplinary across multiple Master’s programs.

A key obstacle to building infectious-disease digital twins is that relevant computational models remain scattered across publications and repositories, described with inconsistent metadata, and difficult to compare, reproduce, or reuse. Yet these models are essential: they provide mechanistic formalisms to represent processes across scales (molecular, cellular, tissue, host, and sometimes population levels), and can be encoded in standard formats such as SBML [4] for quantitative biochemical networks or SBML-qual for qualitative and logical representations of regulatory and signalling systems. Without a structured way to organise these resources, the field loses time reinventing models rather than integrating and extending them.

Beyond long-term research goals, a centralised and well-annotated repository becomes especially critical during emerging outbreaks and health emergencies. In crisis situations, time is a limiting factor : researchers and decision-makers need rapid access to trustworthy models and structured knowledge to understand how a pathogen affects the host, anticipate disease dynamics, and evaluate potential interventions. By providing a single resource where models, metadata, and documentation are gathered in a consistent way. It overcomes the limited reusability of dispersed data, it saves researchers time in model discovery, comparison, and reuse.

To address this gap, InfectioGIT is designed as an open, FAIR[5]-compliant, and community-driven resource that gathers, annotates, and organises computational models related to infectious diseases and host–pathogen interactions. It merges a GitHub repository with a lightweight structured database, creating a unified ecosystem for model discovery, standardised annotation, interoperability, and reuse. By centralising these modelling building blocks and making them easier to curate collaboratively, InfectioGIT aims to provide the scalable digital foundation needed to support infectious-disease digital twins, including One Health–oriented applications such as West Nile virus modelling within the EUR UNITEID framework.

## Host team

This project is supervised by Professor Anna Niarakis and Dr. Virginie Jouffret at the Center for Integrative Biology (CBI) in Toulouse. Professor Niarakis is a specialist in systems biology and leads the Computational Systems Biology team. Dr. Jouffret is an expert in bioinformatics and biostatistics from the BigA platform. Their collaboration at CBI provides the necessary expertise in both biological modeling and data processing for this project.

The CBI is an international multidisciplinary research institute dedicated to understanding biological systems through integrative and multi-scale approaches. It federates several laboratories, including the host laboratory Molecular, Cellular and Developmental Biology (MCD), and addresses major health-related societal challenges by combining experimental biology with quantitative and computational methods.

Within this scientific environment, the host team contributes to the development of integrative computational approaches at the crossroads of mathematics, computer science, and bioinformatics. Their expertise covers the design of reproducible workflows, the structuring of biological knowledge, and the development of tools that make computational resources easier to access, assess, and reuse by the community.

---

# General Objective

The objective of the InfectioGIT project is to collect, organize, standardize, and make accessible computational models describing human infectious diseases, using biomedical ontologies and the BioModels [6] repository as the primary source of mechanistic models.

The project aims to build a centralized, structured, and queryable resource of SBML models relevant to emerging infectious diseases.

---

# Key Steps

## Development of a curated catalogue of infectious diseases

### Selection Criteria

Diseases will be selected based on two complementary approaches:

#### Public health priority selection

Using emergent diseases which refers to diseases whose incidence or geographic range is increasing, that can (re-)appear unexpectedly, and/or that pose a growing outbreak risk due to factors like climate-driven vector expansion, animal–human spillover, travel, and gaps in population immunity, relevant for arboviruses (dengue, chikungunya, West Nile) and other outbreak-prone infections.

- WHO global health priorities  
- National priorities defined by Santé Publique France  
- Pandemic preparedness frameworks  

#### Algorithmic prioritization

Through high-modeling diseases, those with an unusually large number of published computational models and rich datasets, often because they have had major outbreaks, sustained surveillance, and strong research investment (e.g., COVID-19, dengue). Including them is useful for benchmarking pipelines and comparing model types at scale.

A Python-based scoring algorithm is to be developed to rank diseases based on:

- Environmental exposure factors  
- Host diversity  
- Scarcity of existing computational models  
- Urgency of epidemiological reporting  

### Selected Diseases (initial panel)

**National (France) :**

- West Nile virus  
- Dengue [7]  
- Chikungunya [8]  

**International:**

- Mpox [9]  
- Avian influenza (bird flu) and Measles Influenza [10]  
- COVID-19 [11]  

---

## Model Retrieval, Filtering, and Metadata Standardization Pipeline

To ensure the construction of a robust, comprehensive, and reusable repository, the InfectioGIT project will implement a structured pipeline dedicated to the retrieval, filtering, and standardization of computational models and their associated metadata.

### Ontology-Based Data Extraction

To extract relevant data, we will primarily rely on the Disease Ontology (DO) [12], a standardized hierarchical classification of human diseases providing unique identifiers (DOID). This ontology serves as the semantic backbone of the project, enabling precise disease identification and eliminating ambiguity caused by inconsistent naming conventions across databases.

For each selected disease, two complementary ontology identifiers will be collected:

- **DOID identifiers**, extracted from the HumanDO.obo file using a custom Python parsing script  
- **MONDO identifiers** [13], obtained from the MONDO ontology, which integrates multiple disease vocabularies into a unified and interoperable framework  

This dual-ontology strategy is essential to maximize model retrieval coverage. Indeed, BioModels entries are inconsistently annotated: some rely exclusively on DOID, others on MONDO, and some only include textual disease references. By systematically combining both identifiers, we ensure a comprehensive and non-redundant collection of relevant models.

Example of extracted identifier files:

/IDs/

├── dengue_DOID.txt

├── dengue_MONDO.txt

├── chikungunya_DOID.txt

└── chikungunya_MONDO.txt

---

### Definition of Search Keywords

For each disease, a standardized list of search terms will be defined to improve retrieval sensitivity and account for annotation variability within BioModels.

This list will include:

- Synonyms (e.g., alternative disease names)  
- Scientific names of pathogens  
- Alternative spellings  
- Acronyms and abbreviations  

This step is critical because model annotations are not uniformly standardized, and relevant models may only be discoverable through textual queries rather than ontology identifiers.

---

### Retrieval of Models from BioModels

#### Data Access

Models will be retrieved programmatically using:

- The **BioModels REST API**  
- The **Python Bioservices package**, which provides a convenient interface for querying biological databases  

#### Retrieval Pipeline

A dedicated Python script will be developed to automate the retrieval process. The pipeline will:

- Query BioModels using:
  - DOID identifiers  
  - MONDO identifiers  
  - keyword lists  

- Aggregate and deduplicate results obtained from multiple query strategies  

- Filter models based on initial relevance criteria (e.g., presence of associated publication)  

- Download corresponding **SBML files**  

Example of generated output structure:

/models/dengue/

├── BIOMDxxxx.xml

├── metadata_dengue.json

└── model_list.txt

This automated workflow ensures scalability, reproducibility, and consistency in data acquisition.

---

### State of the Art: Disease-Specific Model Landscape

A systematic state-of-the-art analysis will be conducted to refine the selection of approximately **30 to 60 models**. This analysis will characterize the modeling landscape for each disease and guide informed selection.

#### Quantitative Assessment (Model Volume)

The number of available SBML and ODE-based models in BioModels will be quantified for each disease.

This step aims to:

- Identify **high-modeling diseases** (e.g., COVID-19, Dengue)  
- Provide benchmarks for evaluating the retrieval pipeline  
- Assess the balance between well-studied and underrepresented diseases  

#### Redundancy Management

In cases where a large number of models are available for a given disease, a selection strategy will be applied to avoid redundancy while preserving diversity.

Selection will prioritize:

- Scientific relevance  
- Diversity of modeling approaches (e.g., mechanistic, logical, epidemiological)  

---

### Model Filtering and Selection

To ensure the scientific quality of the repository, candidate models will be evaluated through cross-referencing with bibliographic databases such as PubMed and Google Scholar.

#### Selection Criteria

Models will be filtered based on:

- Biological relevance  
- Availability of an associated scientific publication  
- Annotation completeness  
- Reproducibility and clarity of model structure  

#### Citation-Based Filtering

A Python script leveraging **Biopython’s Entrez module** will retrieve publication metadata and apply quantitative filters based on:

- Number of citations  
- Journal quality  
- Publication year  

Priority will be given to widely recognized and validated models (i.e., highly cited publications), ensuring that the repository is built upon robust and credible scientific foundations.

---

### Identification of Knowledge Gaps

An important objective of the project is to identify and document gaps in the current modeling landscape.

Two main types of gaps will be investigated:

- **Scarcity gaps**: diseases with high epidemiological relevance but few available computational models (e.g., Mpox, West Nile virus)  
- **Biological gaps**: absence of models incorporating key factors such as host diversity, environmental exposure, or multi-scale interactions  

These insights will help guide future modeling efforts and highlight areas where additional research is needed.

---

### Metadata Acquisition, Structuring, and Standardization

#### Metadata Sources

Metadata will be collected and integrated from multiple sources:

- BioModels metadata  
- COMBINE archive metadata [14]  
- Disease Ontology (DOID)  
- MONDO ontology  

#### Metadata Standard

A unified metadata schema will be defined to ensure consistency and interoperability across the repository.

This schema will be based on:

- COMBINE Archive standards  
- JSON structured format  
- FAIR data principles  

Example of metadata structure:

```json
{
  "model_id": "BIOMD0000000957",
  "disease": "Dengue",
  "ontology_ids": ["DOID:12365", "MONDO:0001234"],
  "model_type": "SIR",
  "organism": "Human",
  "publication": "...",
  "authors": "...",
  "year": 2020
}
```

This standardized format enables efficient indexing, querying, and interoperability with external bioinformatics tools and databases.

#### Automated Metadata Enrichment

To address inconsistencies and missing annotations in non-curated models, an automated metadata enrichment pipeline will be implemented.

##### Enrichment Process

- **Gap Identification**
  
  A Python script will detect models with missing or incomplete ontology links (DOID/MONDO)

- **Automated Injection**
  
  The script will parse ontology source files to retrieve:
    - Synonyms
    - Parent terms
    - Cross-references
    - Metadata Update

These attributes will be automatically injected into the local metadata.json files

This process ensures that all models reach a consistent level of annotation quality and remain fully searchable within the InfectioGIT ecosystem.

---

# Bibliography

1. Niarakis, A., Laubenbacher, R., An, G., Ilan, Y., Fisher, J., Flobak, Å., Reiche, K., Rodríguez Martínez, M., Geris, L., Ladeira, L., Veschini, L., Blinov, M. L., Messina, F., Fonseca, L. L., Ferreira, S., Montagud, A., Noël, V., Marku, M., Tsirvouli, E., … Glazier, J. A. (2024). Immune digital twins for complex human pathologies: applications, limitations, and challenges. NPJ Systems Biology and Applications, 10(1), 141. https://doi.org/10.1038/s41540-024-00450-5

2. Sinclair, J. R. (2019). Importance of a One Health approach in advancing global health security and the Sustainable Development Goals. Revue Scientifique et Technique (International Office of Epizootics), 38(1), 145–154. https://doi.org/10.20506/rst.38.1.2949

3. Martin, M.-F., & Simonin, Y. (2019). [West Nile virus historical progression in Europe]. Virologie (Montrouge, France), 23(5), 265–270. https://doi.org/10.1684/vir.2019.0787

4. Chaouiya, C., Bérenguier, D., Keating, S. M., Naldi, A., van Iersel, M. P., Rodriguez, N., Dräger, A., Büchel, F., Cokelaer, T., Kowal, B., Wicks, B., Gonçalves, E., Dorier, J., Page, M., Monteiro, P. T., von Kamp, A., Xenarios, I., de Jong, H., Hucka, M., … Helikar, T. (2013). SBML qualitative models: a model representation format and infrastructure to foster interactions between qualitative modelling formalisms and tools. BMC Systems Biology, 7, 135. https://doi.org/10.1186/1752-0509-7-135

5. Wilkinson, M. D., Dumontier, M., Aalbersberg, I. J. J., Appleton, G., Axton, M., Baak, A., Blomberg, N., Boiten, J.-W., da Silva Santos, L. B., Bourne, P. E., Bouwman, J., Brookes, A. J., Clark, T., Crosas, M., Dillo, I., Dumon, O., Edmunds, S., Evelo, C. T., Finkers, R., … Mons, B. (2016). The FAIR Guiding Principles for scientific data management and stewardship. Scientific Data, 3, 160018. https://doi.org/10.1038/sdata.2016.18

6. Le Novère, N., Bornstein, B., Broicher, A., Courtot, M., Donizelli, M., Dharuri, H., Li, L., Sauro, H., Schilstra, M., Shapiro, B., Snoep, J. L., & Hucka, M. (2006). BioModels Database: a free, centralized database of curated, published, quantitative kinetic models of biochemical and cellular systems. Nucleic Acids Research, 34(Database issue), D689-91. https://doi.org/10.1093/nar/gkj092

7. Guzman, M. G., & Harris, E. (2015). Dengue. The Lancet, 385(9966), 453–465. https://doi.org/10.1016/S0140-6736(14)60572-9

8. Maure, C., Khazhidinov, K., Kang, H., Auzenbergs, M., Moyersoen, P., Abbas, K., Santos, G. M. L., Medina, L. M. H., Wartel, T. A., Kim, J. H., Clemens, J., & Sahastrabuddhe, S. (2024). Chikungunya vaccine development, challenges, and pathway toward public health impact. Vaccine, 42(26), 126483. https://doi.org/10.1016/j.vaccine.2024.126483

9. Hou, W., Wu, N., Liu, Y., Tang, Y., Quan, Q., Luo, Y., & Jin, C. (2025). Mpox: Global epidemic situation and countermeasures. Virulence, 16(1), 2457958. https://doi.org/10.1080/21505594.2025.2457958

10. Bi, Y., Yang, J., Wang, L., Ran, L., & Gao, G. F. (2024). Ecology and evolution of avian influenza viruses. Current Biology, 34(15), R716–R721. https://doi.org/10.1016/j.cub.2024.05.053

11. Ochani, R., Asad, A., Yasmin, F., Shaikh, S., Khalid, H., Batra, S., Sohail, M. R., Mahmood, S. F., Ochani, R., Hussham Arshad, M., Kumar, A., & Surani, S. (2021). COVID-19 pandemic: from origins to outcomes. A comprehensive review of viral pathogenesis, clinical manifestations, diagnostic evaluation, and management. Le Infezioni in Medicina : Rivista Periodica Di Eziologia, Epidemiologia, Diagnostica, Clinica e Terapia Delle Patologie Infettive, 29(1), 20–36.

12. Schriml, L. M., Mitraka, E., Munro, J., Tauber, B., Schor, M., Nickle, L., Felix, V., Jeng, L., Bearer, C., Lichenstein, R., Bisordi, K., Campion, N., Hyman, B., Kurland, D., Oates, C. P., Kibbey, S., Sreekumar, P., Le, C., Giglio, M., & Greene, C. (2019). Human Disease Ontology 2018 update: classification, content and workflow expansion. Nucleic Acids Research, 47(D1), D955–D962. https://doi.org/10.1093/nar/gky1032

13. Vasilevsky, N. A., Toro, S., Matentzoglu, N., Flack, J. E., Mullen, K. R., Hegde, H., Gehrke, S., Whetzel, P. L., Shwetar, Y., Harris, N. L., Ngu, M. S., Alyea, G. L., Kane, M. S., Roncaglia, P., Sid, E., Thaxton, C. L., Wood, V., Abraham, R. S., Achatz, M. I., … Haendel, M. A. (2025). Mondo: integrating disease terminology across communities. Genetics. https://doi.org/10.1093/genetics/iyaf215

14. Gennari, J. H., König, M., Misirli, G., Neal, M. L., Nickerson, D. P., & Waltemath, D. (2021). OMEX metadata specification (version 1.2). Journal of Integrative Bioinformatics, 18(3). https://doi.org/10.1515/jib-2021-0020

---
