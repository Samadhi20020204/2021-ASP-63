Differential Gene Expression Analysis Between Cancer and Normal Tissue Samples Using R

1. Title and Abstract

Project Title

Differential Gene Expression Analysis Between Cancer and Normal Tissue Samples Using R

2. Introduction

2.1 Background

Gene expression is the process through which information stored in genes is used to produce functional products such as proteins. Gene expression levels can vary between different cell types and biological conditions.

Changes in gene expression may be associated with several diseases, including cancer. Therefore, analyzing gene expression data can help researchers understand biological differences between healthy and diseased tissues.

2.2 Research Problem

Cancer is associated with changes in normal cellular processes, including gene expression. Identifying genes that show different expression levels between cancer and normal tissue samples is an important step in understanding these changes.

However, gene expression datasets may contain a large number of genes and samples. Analyzing these datasets manually can be difficult. Bioinformatics tools and statistical methods can help identify meaningful expression patterns.

2.3 Research Aim

The main aim of this project is to identify and analyze differentially expressed genes between cancer and normal tissue samples using bioinformatics techniques in R.

2.4 Research Objectives

To obtain a publicly available gene expression dataset containing cancer and normal tissue samples.

To inspect and preprocess the selected dataset using R.

To explore gene expression patterns through data visualization.

To identify genes with significant expression differences between cancer and normal tissue samples.

To interpret and present the findings using appropriate statistical methods and visualizations.

2.5 Research Question

Which genes show significant differences in expression between cancer and normal tissue samples?

2.6 Hypothesis

Null Hypothesis (H0): There are no significant differences in gene expression between cancer and normal tissue samples.

Alternative Hypothesis (H1): There are significant differences in gene expression between cancer and normal tissue samples.

3. Data

3.1 Dataset Source

The gene expression dataset used in this project was obtained from the Gene Expression Omnibus (GEO), a public repository maintained by the National Center for Biotechnology Information (NCBI).

The selected dataset is GSE15852.

Dataset link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15852

3.2 Dataset Description

GSE15852 is a human breast cancer gene expression dataset containing 43 breast tumor samples and their 43 paired normal tissue samples.

The dataset contains a total of 86 samples from 43 patients, with one cancer sample and one corresponding normal tissue sample from each patient.

Gene expression measurements were generated using the Affymetrix Human Genome U133A microarray platform (GPL96).

The dataset contains 22,283 probe sets and 86 samples, resulting in an expression matrix of 22,283 × 86.

This dataset will be used to compare gene expression levels between breast cancer and normal tissue samples.

3.3 Variables

The main variables considered in this project include:

- **Probe Set ID:** Identifier of each measured probe set on the microarray platform.
- **Gene Expression Values:** Expression measurements obtained from the microarray experiment.
- **Sample ID:** Unique identifier assigned to each sample.
- **Sample Group:** Indicates whether a sample belongs to the Cancer or Normal group.
- **Patient ID:** Identifies the patient associated with each sample.
- **Patient Pair:** Identifies the corresponding cancer and normal samples collected from the same patient.

During data inspection, the expression matrix was confirmed to contain 22,283 probe sets and 86 samples. The sample metadata was verified to contain 43 cancer samples and 43 normal samples from 43 paired patients. The sample order in the expression matrix was also confirmed to match the sample metadata.

3.4 Ethical and Licensing Considerations

This project uses publicly available gene expression data from the NCBI GEO database. The original dataset and its associated publication will be acknowledged in the references.

No new biological samples or personal information will be collected for this project. The dataset provider's usage conditions will be followed.

3.5 Data Acquisition and Initial Preparation

The GSE15852 dataset was downloaded from the NCBI Gene Expression Omnibus (GEO) database using the GEOquery package in R.

The downloaded dataset was imported into R as an ExpressionSet object. The gene expression measurements and sample metadata were extracted separately for further analysis.

The sample metadata was processed to identify the sample groups and corresponding patient IDs. The dataset was verified to contain 43 cancer samples and 43 paired normal samples.

The expression data and processed sample metadata were saved in RDS and CSV formats for subsequent analysis.
