# PGPR Gene Detector: A Protein Sequence–Based Computational Tool for PGPR Gene Identification
---
## Introduction

Plant Growth-Promoting Rhizobacteria (PGPR) are beneficial soil microorganisms that significantly contribute to plant health and productivity. These bacteria facilitate plant growth through various direct and indirect mechanisms, making them important targets for agricultural biotechnology and microbial research.

Traditional approaches for PGPR gene identification, such as manual consensus sequence generation using tools like BioEdit, become inefficient when applied to large datasets or long gene sequences. Additionally, nucleotide-level comparisons often fail to accurately reflect functional similarity.  

This project introduces a protein-level computational framework that improves scalability, accuracy, and biological relevance by leveraging sequence alignment–based analysis.

---

## Methodology

### Data Collection

PGPR-related gene and protein sequences were collected from the National Center for Biotechnology Information (NCBI) database, along with curated datasets contributed by previous academic research efforts. The initial dataset consisted of approximately 1000 sequences.

### Data Preprocessing

To ensure data quality and reliability, sequences were filtered to remove:
- Incomplete records  
- Discontinued or outdated entries  
- Low-quality or redundant sequences  

After preprocessing, a final dataset of **265 validated PGPR protein sequences** was retained for analysis.

### Protein-Level Analysis Rationale

The analysis was conducted at the protein level due to the following advantages:
- Protein sequences are shorter and computationally efficient
- Functional similarity is better conserved at the protein level
- Reduced noise compared to nucleotide-level variation

---

## System Workflow

1.	User submits a protein sequence  
2. The sequence is aligned against the PGPR dataset using pairwise alignment  
3. Similarity scores are calculated for all reference sequences  
4. The top 10 most similar PGPR sequences are identified  
5. Classification is performed based on a similarity threshold (default ≥ 70%)  
6. A consensus sequence is generated from top matches  
---

## Features

### PGPR Sequence Analyzer
- Accepts protein sequence input  
- Performs pairwise local alignment using Biopython  
- Displays top 10 similar PGPR sequences with similarity scores  
- Classifies sequences as PGPR or Non-PGPR  
- Generates a consensus sequence for biological interpretation  

### DNA to Protein Converter
If you have DNA sequence first convert it then Analysis.
- Validates DNA sequences (A, T, G, C only)  
- Ensures sequence length is divisible by three  
- Translates valid DNA sequences into protein sequences  
- Provides informative error messages for invalid input  

### Add New Gene to Dataset
- Allows users to submit new PGPR gene entries 
- Submissions are stored temporarily   
- Administrative approval is required before inclusion in the main dataset  

### Dataset Download
- Provides access to the validated PGPR dataset in CSV format  
- Enables offline analysis and reproducibility  

---

## Implementation Details

- **Programming Language:** Python  
- **Web Framework:** Flask  
- **Bioinformatics Library:** Biopython  
- **Alignment Strategy:** Pairwise Local Alignment  
- **Data Storage:** CSV-based datasets  

---

## Deployment

The application is deployed using **Render**, with continuous integration enabled through GitHub. Any updates pushed to the repository are automatically reflected in the deployed application, ensuring ease of maintenance and scalability.

---

## Results and Validation

- Known PGPR protein sequences (e.g., *Pseudomonas protegens*) were correctly classified as PGPR  
- Non-PGPR protein sequences (e.g., *Escherichia coli* cell division proteins) were correctly classified as non-PGPR  
- The tool provides ranked similarity scores, final classification, and consensus sequence output  

---

## Limitations

- Classification is based on a fixed similarity threshold  
- Performance depends on dataset size and diversity  
- Currently supports single-sequence input only  

---

## Future Scope

- Integration of machine learning–based classification models  
- Batch FASTA file upload support    
- Expansion of PGPR gene families  

---

## Acknowledgements

This project was developed as part of an academic bioinformatics initiative under the guidance ofProf. Ashutosh Mani. The work was carried out by Himanshu Singhal and Samarth Beck, who served as primary and main contributors with collaborative support from the project team.
---

## License

This software is released for **academic and non-commercial research purposes only**.
Commercial use requires prior permission from the authors.
