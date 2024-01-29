# Histone marks identify novel transcription factors that parse CAR-T cell subset-of-origin, clinical potential and expansion

## Content

- [Overview](##overview)
- [Repository Aims](##repository-aims)
- [System Requirements](##system-requirements)
- [What is CUT&RUN](##what-is-cut&run)
- [Data Processing and Analysis Pipeline](##workflow)
- [Raw Data Processing](##raw-cut&run-data-processing)
- [Quality Check and Basic Analysis](##quality-check-and-basic-analysis)
- [Analysis in the Manuscript and Figures Generation](#analysis-in-the-manuscript-and-figures-generation)
- [Contact](#contact)
- [License](./LICENSE)

## Overview

Chimeric antigen receptor-modified T cell (CAR-T) immunotherapy has revolutionized the treatment of blood cancers. Parsing the genetic underpinnings of T cell precursor quality and subsequent CAR-T efficacy is challenging. RNA-seq informs infused CAR-T state, but the nature of dynamic transcription during activation hinders identification of transiently or minimally expressed genes, such as transcription factors, and over-emphasizes effector and metabolism genes. We investigated whether analyses of transcriptionally repressive and permissive histone methylation marks reveal associations with CAR-T potential beyond what is seen by transcriptomic analysis. We assessed human CD8+ T cell naïve, central and effector memory subsets that form the substrate of CAR-T cell products, and CAR-T cells derived from these subsets. We extended these observations into the clinic, by examining CAR-T products from a clinical trial of lymphoma patients (NCT01865617). We report that histone marks provide a rich dataset for identification of genes not apparent by conventional transcriptomics. Histone marks improved identification of T cell subsets, CAR-T manufactured from these subsets, and CAR-T manufactured from central memory cells from healthy donors and patients. Using this discriminative approach, we controlled for clinical factors and identified a factor, KLF7, associated with CAR-T cell expansion in patients. Epigenomic methods are an orthogonal, robust and wide-reaching approach for the assessment of T cell immunotherapeutic quality.

## Repository Aims

The primary aim of this GitHub repository is to facilitate the reproduction of the results presented in our manuscript, "Histone marks identify novel transcription factors that parse CAR-T cell subset-of-origin", clinical potential and expansion. We provide reproducible codes and description for CUT&RUN and RNA-seq data processing, quality check, basic analysis and customized analysis used in the manuscripts, as well as the scripts used for figures generation as outlined in the [Content](#content).

For any questions or further clarifications, users are invited to open an issue within the repository or [contact the authors](##contact) directly.


## What is CUT&RUN

**Cleavage Under Targets and Release Using Nuclease**, **CUT&RUN** for short, is an antibody-targeted chromatin profiling method to measure the histone modification enrichment or transcription factor binding. This is a more advanced technology for epigenomic landscape profiling compared to the tradditional ChIP-seq technology and known for its easy implementation and low cost.  The procedure is carried out in situ where micrococcal nuclease tethered to protein A binds to an antibody of choice and cuts immediately adjacent DNA, releasing DNA-bound to the antibody target. Therefore, CUT&RUN produces precise transcription factor or histone modification profiles while avoiding crosslinking and solubilization issues. Extremely low backgrounds make profiling possible with typically one-tenth of the sequencing depth required for ChIP-seq and permit profiling using low cell numbers (i.e., a few hundred cells) without losing quality.

<img src="./figures/cutrun.png" alt="CUT&RUN (Image source: Skene & Henikoff, ELife 2017 Figure 1A)" width="600px">


## Workflow

**T Cells and CAR-T Cells CUT&RUN Data Processing and Analysis Pipeline**:

<img src="./figures/bioinformatics_workflow.png" alt="Workflow" width="600px">

A more general CUT&RUN or CUT&Tag data processing and analysis tutorial can be found at [Procotols.io: CUT&Tag Data Processing and Analysis Tutorial](https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1) where demo data, demo data analysis codes, and results intepretation are introduced.  

## System Requirements

To run the scripts and analyses contained in this repository, users will need a computer with an operating system capable of running R and shell scripts, such as Windows, macOS, or Linux. The R environment should be set up with a version of R (we recommend R version 4.0 or higher) and Jupyter Notebook, if intending to run R notebooks. For R package installation, users can execute the following command in the R console to install the necessary packages: such as ```install.packages(c("ggplot2"))```, replacing "ggplot2" with the names of the required packages. All the required packages used in the scripts or notebooks are listed at the beginning of each script. Additionally, for running shell scripts, a Unix-like environment is required for Windows users, which can be set up using software such as Git BASH or Windows Subsystem for Linux (WSL). Ensure that all dependencies and required libraries are installed and properly configured before running the scripts to avoid execution errors. All the required software or modules used in the shell scripts are listed at the beginning of each script.


## Raw Data Processing:

- [Processing](./Processing/)

**RNA**

0. [Quality check before alignment](./Processing/RNA_step0_preAlignment_qualitycheck_fastqc.sh)

1. [Alignment](./Processing/RNA_step1_RSEMalignment.sh)

**CUT&RUN**

0. [Preparation before alignment](./Processing/CUTRUN_step0_collectDataRename.sh)

1. [Alignment](./Processing/CUTRUN_step1_alignment.sh)

2. [Peak Calling](./Processing/CUTRUN_step2.1_peakCalling_controlPart.sh)

3. [Post-alignment summary](./Processing/CUTRUN_step3.1_analysis_AlignmentSummary.Rmd)

## Quality Check and Basic Analysis: 

- [Vignettes](./Vignettes/)

**RNA**

1. [Quality Check](./Vignettes/RNAseq_QC.Rmd)

2. [Differentially Expressed Genes Detection](./Vignettes/RNAseq_DE_limma.Rmd)

**CUT&RUN**

1. [Quality Check](./Vignettes/CUTRUN_QC.Rmd)

2. [Differential Enriched Histone Detection](./Vignettes/CUTRUN_DifferentialDetection.Rmd)

3. [Peak Calling and Comparison Analysis](./Vignettes/CUTRUN_PeakCalling.Rmd)

## Analysis in the Manuscript and Figures Generation: 

- [Notebooks](./Notebooks/)

**Quality Check**

    - [Paper_Figures_SupplementaryQC.ipynb](./Notebooks/Paper_Figures_SupplementaryQC.ipynb)

    - [Paper_Figures_SupplementaryQC_RNA.ipynb](./Notebooks/Paper_Figures_SupplementaryQC_RNA.ipynb)

**Multi-omics Joint Analysis**
    
    - [Paper_Figures_CUTRUN_RNA_JointAnalysis.ipynb](./Notebooks/Paper_Figures_CUTRUN_RNA_JointAnalysis.ipynb)

    - [Paper_Figures_CUTRUN_RNA_Comparison.ipynb](./Notebooks/Paper_Figures_CUTRUN_RNA_Comparison.ipynb)


**Differential Analysis**

    - [Paper_Figures_DifferentialEnrichedHistoneAssociatedGenes.ipynb](./Notebooks/Paper_Figures_DifferentialEnrichedHistoneAssociatedGenes.ipynb)

    - [Paper_Figures_DifferentialHistoneEnrichmentAnalysis.ipynb](./Notebooks/Paper_Figures_DifferentialHistoneEnrichmentAnalysis.ipynb)

**Bivalent Genes**

    - [Paper_Figures_HistoneRNA_PCA_Bivalency.ipynb](./Notebooks/Paper_Figures_HistoneRNA_PCA_Bivalency.ipynb)

    - [Paper_Figures_BivalencyScatterPlot_ClinicalDifferentialAnalysis.ipynb](./Notebooks/Paper_Figures_BivalencyScatterPlot_ClinicalDifferentialAnalysis.ipynb)


**Patient Clinical Associated Analysis**

    - [Paper_Figures_PatientClinicalOutcomesDifferentialAnalysis.ipynb](./Notebooks/Paper_Figures_PatientClinicalOutcomesDifferentialAnalysis.ipynb)

    - [Paper_Figures_PatientSamples_ExpansionAnalysis.ipynb](./Notebooks/Paper_Figures_PatientSamples_ExpansionAnalysis.ipynb)

    - [Paper_Figures_CM-derivedCART_vs_CM-healthyDonor_CUTRUN.ipynb](./Notebooks/Paper_Figures_CM-derivedCART_vs_CM-healthyDonor_CUTRUN.ipynb)

    - [Paper_Figures_CM-derivedCART_vs_CM-healthyDonor_RNA.ipynb](./Notebooks/Paper_Figures_CM-derivedCART_vs_CM-healthyDonor_RNA.ipynb)

## Contact

[Ye Zheng](https://yezhengstat.github.io/)

yzheng23@fredhutch.org

## License
This code repository is covered under GNU General Public License (GPL) v3.0.