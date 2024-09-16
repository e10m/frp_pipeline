# Fragment Recruitment Plot Pipeline

**Author:** Dien Mach  
**UNIX (GCP) Username:** ethanmach1998  
**GitHub Repository:** [https://github.com/e10m/frp_pipeline](https://github.com/e10m/frp_pipeline)

---

## Overview

The Fragment Recruitment Plot Pipeline is a bioinformatics project that involves gathering metagenomic reads, aligning them to reference genomes using two different alignment algorithms, and generating Fragment Recruitment Plots (FRPs) to visualize the alignment results. This project demonstrates proficiency in bioinformatics tools, workflow automation, and data visualization, providing insights into microbial community composition.

> **Note:** This README provides a comprehensive overview of the project's functionality, design, and implementation to showcase my work to potential employers.

---

## Table of Contents

- [Data Gathering](#data-gathering)
- [Technologies Used](#technologies-used)
- [Project Motivation](#project-motivation)
- [Implementation Details](#implementation-details)
  - [Alignment to Reference Genomes](#alignment-to-reference-genomes)
  - [Fragment Recruitment Plots (FRPs)](#fragment-recruitment-plots-frps)
  - [Workflow Automation with Snakemake](#workflow-automation-with-snakemake)
- [Analysis Highlights](#analysis-highlights)
- [Challenges and Solutions](#challenges-and-solutions)
- [Learning Outcomes](#learning-outcomes)
- [Future Enhancements](#future-enhancements)
- [References](#references)
- [Contact Information](#contact-information)

---

## Data Gathering

- **Metagenomic Sample**: Used the Human Microbiome Project (HMP) mock community reads obtained using Roche 454 sequencing technology.
- **Files**: Downloaded FASTQ and FASTA files for 16S rRNA reads.
- **Reference Genomes**: Selected 10 bacterial reference genomes from the NCBI database to align reads against.
- **Storage**: Files related to this project are stored on Google Cloud in the `AS41073481-dmach1` project within the `metagenomics` instance.

### Reference Genomes Table

| Organism Name & Repository Number     | Taxon ID | Contig/Scaffold Count | GenBank Gene Count |
|---------------------------------------|----------|-----------------------|--------------------|
| Acinetobacter baumannii ATCC 17978    | 470      | 2/2                   | 3,839              |
| Staphylococcus epidermidis ATCC 12228 | 1282     | 2/2                   | 2,354              |
| Streptococcus pneumoniae ATCC BAA-334 | 170187   | 1/1                   | 2,292              |
| Deinococcus radiodurans DSM 20539     | 243230   | 4/4                   | 3,305              |
| Lactobacillus gasseri DSM 20243       | 324831   | 26                    | 1,786              |
| Actinomyces odontolyticus ATCC 17982  | 411466   | 4/2                   | 2,214              |
| Methanobrevibacter smithii ATCC 35061 | 420247   | 1/1                   | 1,837              |
| Bacteroides vulgatus ATCC 8482        | 435590   | 1/1                   | 4,183              |
| Enterococcus faecalis ATCC 47077      | 474186   | 75/NA                 | 2,584              |
| Bacillus cereus ATCC 10987            | 2026187  | 239/222               | 5,807              |

*Table 1: Reference genomes downloaded from the NCBI website.*

---

## Technologies Used

- **Programming Languages:**
  - **Python**
    - `matplotlib`: Data visualization library.
    - `plotly`: Interactive graphing library.
    - `re`: Regular expressions library for parsing.
- **Bioinformatics Tools:**
  - **Bowtie2**: Tool for aligning sequencing reads.
  - **BWA-MEM**: Software package for mapping sequences against large reference genomes.
  - **SAMtools**: Suite for interacting with alignment data.
- **Workflow Management:**
  - **Snakemake**: Workflow management system for automating pipelines.
- **Package Managers:**
  - **Anaconda**: Python distribution for data science.
- **Scripting and Automation:**
  - **Bash scripting** for executing Snakemake workflows.

---

## Project Motivation

The goal of this project was to:

- Analyze metagenomic data to determine the composition of microbial communities.
- Compare the performance of different alignment algorithms (Bowtie2 and BWA-MEM).
- Develop an automated pipeline to streamline the data processing workflow.
- Gain proficiency in bioinformatics tools and data visualization techniques.

---

## Implementation Details

### Alignment to Reference Genomes

- **Alignment Algorithms Used:**
  - **Bowtie2**
  - **BWA-MEM**

- **Standardizing Reference Genomes:**
  - Developed a Python script, `contig_combiner.py`, to concatenate contigs/scaffolds in the FASTA files.
  - Artificially linked contigs using a nucleotide sequence (`XXXXX`) to create a continuous reference for alignment.

- **Process:**
  - Indexed the reference genomes using both Bowtie2 and BWA-MEM indexing tools.
  - Aligned the metagenomic reads to each reference genome using both algorithms.
  - Generated SAM files containing alignment data for further analysis.

### Fragment Recruitment Plots (FRPs)

- **Purpose:** Visualize the alignment of metagenomic reads to reference genomes to assess the presence and abundance of organisms in the sample.

- **Script Developed:** `sam_to_frp.py`
  - Parses SAM files generated from alignments.
  - Calculates read coverage across the reference genome.
  - Generates FRPs using **Matplotlib** and **Plotly** for interactive visualization.

- **Example FRPs:**

  ![FRPs for Taxon ID 243230](./images/243230_frps.png)

  *Figure 1: FRPs for Taxon ID 243230 (*Deinococcus radiodurans* DSM 20539) aligned using Bowtie2 (top) and BWA-MEM (bottom).*

### Workflow Automation with Snakemake

- **Automation:**
  - Utilized **Snakemake** to automate the alignment and FRP generation processes.
  - Wrote a Bash script, `main.sh`, containing all Snakemake commands to generate 20 FRPs efficiently.

- **Workflow Diagram:**

  ![Snakemake Workflow Diagram](./images/frp_pipeline.png)

  *Figure 2: Snakemake workflow diagram illustrating the Directed Acyclic Graph (DAG) of the pipeline.*

- **Benefits:**
  - Streamlined the workflow, reducing manual intervention.
  - Improved reproducibility and scalability of the analysis.
  - Reduced processing time by approximately **60%** compared to manual execution.

---

## Analysis Highlights

- **Most Represented Genome:**
  - **Taxon ID 243230** (*Deinococcus radiodurans* DSM 20539) with a **35.24% overall alignment rate**.

- **Algorithm Comparison:**
  - **Bowtie2:**
    - Faster runtimes.
    - Provided more descriptive alignment statistics.
  - **BWA-MEM:**
    - Slower but may offer different alignment sensitivity.

- **Findings:**
  - The high representation of *Deinococcus radiodurans* suggests its abundance in the mock community sample.
  - Differences in alignment results highlight the importance of selecting appropriate algorithms for metagenomic analysis.

---

## Challenges and Solutions

- **Workflow Complexity:**
  - **Challenge:** Managing multiple alignment tasks and plot generations manually was inefficient.
  - **Solution:** Implemented **Snakemake** to automate the pipeline, significantly reducing processing time.

- **Standardizing Reference Genomes:**
  - **Challenge:** Variability in contig/scaffold formats across different reference genomes.
  - **Solution:** Developed `contig_combiner.py` to standardize references, ensuring consistent alignment.

- **Data Handling:**
  - **Challenge:** Handling large genomic datasets and ensuring computational efficiency.
  - **Solution:** Optimized scripts and leveraged efficient data structures for parsing and processing.

---

## Learning Outcomes

- **Technical Proficiency:**
  - Gained experience with bioinformatics tools like Bowtie2, BWA-MEM, SAMtools, and Snakemake.
  - Improved Python scripting skills, particularly in data parsing and visualization.

- **Workflow Automation:**
  - Learned to automate complex bioinformatics workflows, enhancing reproducibility and efficiency.

- **Data Analysis and Interpretation:**
  - Developed skills in interpreting alignment results and generating meaningful visualizations.

- **Problem-Solving:**
  - Overcame challenges related to data standardization and visualization.

---

## Future Enhancements

- **Expand Reference Genomes:**
  - Include more diverse genomes to broaden the analysis scope.

- **Algorithm Optimization:**
  - Explore additional alignment algorithms or parameters to improve accuracy.

- **Enhanced Visualization:**
  - Develop interactive web-based visualizations for easier data exploration.

- **Statistical Analysis:**
  - Implement statistical methods to quantify differences between alignment algorithms.

- **Data Visualization with Seaborn or R:**
  - Experiment with different visualization techniques using Seaborn or R to better represent genome positions on the X-axis of each graph.

---

## References

1. **Bowtie2 Documentation:** [http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
2. **BWA-MEM Documentation:** [http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/)
3. **SAMtools Documentation:** [http://www.htslib.org/doc/samtools.html](http://www.htslib.org/doc/samtools.html)
4. **Snakemake Documentation:** [https://snakemake.readthedocs.io/en/stable/](https://snakemake.readthedocs.io/en/stable/)
5. **Matplotlib Documentation:** [https://matplotlib.org/](https://matplotlib.org/)
6. **Plotly Documentation:** [https://plotly.com/python/](https://plotly.com/python/)
7. **NCBI Reference Genomes:** [https://www.ncbi.nlm.nih.gov/refseq/](https://www.ncbi.nlm.nih.gov/refseq/)
8. **Anaconda Distribution:** [https://www.anaconda.com/](https://www.anaconda.com/)
9. **Bioconda:** [https://bioconda.github.io/](https://bioconda.github.io/)
10. **Human Microbiome Project (HMP):** [https://www.hmpdacc.org/](https://www.hmpdacc.org/)

---

## Contact Information

For more information about the project or to discuss potential collaborations, please contact:

**Dien Mach**  
Email: [dienethanmach@gmail.com](mailto:dienethanmach@gmail.com)  
GitHub: [e10m](https://github.com/e10m)

---

**Note to Employers:** This project demonstrates my ability to develop and automate bioinformatics workflows, handle large genomic datasets, and visualize complex data. I welcome the opportunity to discuss this project further and answer any questions you may have.
