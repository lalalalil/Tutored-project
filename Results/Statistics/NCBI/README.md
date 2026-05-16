# NCBI Analysis Outputs

## Files

### 01_file_extensions.png
Bar chart showing the top 15 file types found in NCBI article downloads across disease categories (chikungunya, COVID-19, dengue, HIV, influenza, mpox, tuberculosis, West Nile virus).

Files are classified as either DATA (computational files like XML, JSON, CSV) or METADATA (PDFs, images, documents).

Key finding: NCBI downloads are dominated by metadata files. Computational models are rare.

### 02_pmc6620235_sbml_organisms.csv
Table of SBML (Systems Biology Markup Language) organisms/models found in article PMC6620235.

Columns:
- No.: Sequential identifier
- Organism: Model organism or bacterial strain name
- Type: Classification as Tuberculosis, Influenza, or Other

The article contains 20 total SBML files representing 10 unique organisms:
- 1 Tuberculosis model
- 1 Influenza model
- 8 other organisms (Staphylococcus aureus, Klebsiella pneumoniae, Pseudomonas aeruginosa, Haemophilus influenzae, etc.)

This demonstrates that a single article can contain computational models for multiple organisms, often across different disease categories.
