# Lung Cancer Clinical and Gene Expression Profile Data Analysis
This repository analyzes lung cancer patients' clinical and gene expression profile data. This data was the result of a study to identify NSCLC tumor types. In our analysis, we aim to find out which genes are expressed differently in male and female patients. 

# Files:
-`config.yaml:` This file includes the dataset's path, which specifies where the relevant data is kept.  
-`lung_data.ipynb:` Primary Jupyter Notebook for lung cancer analysis. Users may view the code, findings, and visualizations here.  
-`load_df.py:` This Python script loads the dataset for interactive visualization.  
-`test_dash.py:` This file contains code for an interactive dashboard. Running this script enables users to interact with and observe the analysis results. 

# Installation 
To run this repository locally, follow these steps:
1. At first, clone this repository using this command:   
 `git clone https://github.com/eshita53/Programming-2-Course-Assignment-Submission/tree/main/Final_assignment`
2. Then go to the folder which contains the repository contents.
   `cd Programming-2-Course-Assignment-Submission`
3. The next step is to configure your environment. The conda package manager is used in this tutorial. Please ensure that conda or Miniconda are properly installed and configured on your system.
   Run the following command inside the folder to create the environment.
   `conda env create -f environment.yml`
4. Use the following command to activate the environment:
   `conda activate programming-2`
5. From inside the `programming-2` environment, launch the Jupyter notebook server:
   `jupyter notebook`


# Configuaration File
 
The paths to the datasets are contained in a `config.yaml` file. Ensure the dataset locations match the ones in the YAML file before running. The `config` file includes the following:

- `gene_sample_data.ods`: Information regarding the sample of lung cancer patients.
- `gene_series_information.ods`: Information about the experimental setup of sample collections.
- `microarray_data.txt`: Comprehensive information about the prob set and its expression level.
- `Lung3.metadata.xls`: Clinical metadata of 89 patients.

# Run
Start the Jupyter Notebook in the newly created `programming-2` virtual environment and run the `lung_data.ipynb` file.  
To run the interactive dashboard, use `python3 test_dash.py`.


