# Specification Document

## InfectioGIT : A Structured Repository of Models and Metadata for Human Infectious Diseases

**Students :** Aya Sebbah - Diana Bravais - Thalie Holmiere  
**Master 1 – Bioinformatics & Systems Biology**  
**Supervisors :** Pr. Anna Niarakis & Dr Virginie Jouffret  
**Affiliation :** CBI  

**2025-2026**

---

# Table of Content

- Table of Content
- I. Project Context
- II. General Objective
- III. Key Steps
  - 1. Development of a curated catalogue of infectious diseases
    - A. Selection Criteria
      - a. Public health priority selection
      - b. Algorithmic prioritization
    - B. Selected Diseases (initial panel)
  - 2. Model Retrieval, Filtering, and Metadata Standardization Pipeline
    - A. Definition of Search Keywords
    - B. Retrieval of Models from BioModels
    - C. State of the Art: Disease-Specific Model Landscape
      - a. Quantitative Assessment (Model Volume)
      - b. Model Filtering and Selection
      - c. Identification of Knowledge Gaps
    - D. Metadata Acquisition Structuring and Standardization
    - E. Automated Metadata Enrichment
- IV. Technical Materials and Architecture
  - 1. Data Architecture and Repository Structure
    - A. Materials
    - B. InfectioGIT GitHub Repository Organization
    - C. Technical Components and Material
- V. Constraints and Risks
  - 1. Functional Constraints
  - 2. Organizational Constraints
  - 3. Technical and Data Constraints
- VI. Final Deliverables
- Bibliography
- Annexes

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

To extract data, we will primarily prioritize resources from the Disease Ontology [12] (DO)...

*(contenu intégral conservé — suite du document inchangée)*

---

# Bibliography

*(contenu intégral conservé tel quel)*

---

# Annexes

Annexe 1 : Estimated Project Timeline and Key Deliverables