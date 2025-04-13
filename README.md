##  Code for the paper

The entry point for the analysis code is located in the ```paper``` directory. Below is an explanation of each Python script:

### Data Preprocessing
* ```preprocess_proteome.py``` Preprocesses proteomic and RNA expression data.
* ```preprocess_gwas.py``` Preprocesses GWAS summary statistics.

### Tool Execution
* ```run/run_rez.py``` Run REZ tool.
* ```run/run_ldsc.py``` Run S-LDSC tool.
* ```run/run_magma.py``` Run MAGMA tool.
* ```run/run_dese.py``` Run DESE tool.
* ```run/run_ecs.py``` Run ECS tool.

### Analysis and Visualization

* ```analyze_ome_data.py``` Analyzes and visualizes the correlation of the proteome and transcriptome.
* ```analyze_assoc_gene.py``` Analyzes and visualizes results of disease-associated genes.
* ```analyze_assoc_tissue.py``` Analyzes and visualizes results of disease-associated tissues.

### Supporting Tools and Parameters
* ```para.py``` Stores parameters and paths to resource files.
* ```util.py``` Utility functions.