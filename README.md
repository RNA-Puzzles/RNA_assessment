# Basic RNA 3D structure comparison metrics used in RNA-Puzzles assessment. 
---

This git includes:  

* an RNA structure normalization tool used in RNA-Puzzles.
* the RNA 3D structure comparison metrics used in RNA-Puzzles assessment, including RMSD (all-atom), P value, Interaction Network Fidelity and Deformation Index. 

## Installation
To install the package:    
`git clone https://github.com/RNA-Puzzles/RNA_assessment.git`    
`cd RNA_assessment`    
`python setup.py install`    


## Dependency
This package depends on `biopython`, `rna-tools`   
To install:   
`pip install biopython`，  
`pip install rna-tools`,  
OpenStructure for lDDT calculations.

If need to calculate the Interaction Network Fidelity, it needs to call [`MC-annotate`](https://major.iric.ca/MajorLabEn/MC-Tools.html).    
Please download the binary execution from the website and coordinate the directory for it at the top line `MCAnnotate_bin=` of the mcannotate.py script. 

If need to calculate the Local Distance Difference Test (lDDT), it needs to call [OpenStructure](https://openstructure.org/download).    
The lddt folder used in the current script is in [Google Driver](https://drive.google.com/drive/folders/1ZuugpvBi90LG9nW3RZ9BIbxOMDUqxfmw?usp=sharing), which needs to be fetched and extracted in the same directory as example.ipynb


If need to calculate the Atomic Rotationally Equivariant Scorer (ARES), it needs to call [ARES](https://www.science.org/doi/10.1126/science.abe5650)

If need to calculate the Clash score, it needs to call [MolProbity](http://molprobity.biochem.duke.edu/)
        
         
# Repository Structure
This repository contains the codes necessary to replicate all figures from Puzzles Round V. It includes:

- `The figures_reproducibility.ipynb notebook` : this file can generate all the figures found in the Puzzles, with scores holding the original data for the figures. For details on how to generate data using various metrics, please refer to https://github.com/RNA-Puzzles/RNA_assessment.
- `data`: This folder contains the input files necessary for structure normalization.
- `example`: This folder contains the input and output files for example.ipynb/example.py.
- `RNA_normalizer or RNA_normalizer.egg-info or build`: This folder contains the packaged scripts imported in example.ipynb/example.py
- `example.py or example.ipynb`: The folder contains example.py or example.ipynb, which are specific runnable examples of Python scripts or Jupyter notebooks.

# how to use
A detailed introduction can be found in the example [notebook](https://github.com/RNA-Puzzles/RNA_assessment/blob/master/example.ipynb) or the example [script](https://github.com/RNA-Puzzles/RNA_assessment/blob/master/example/example.py). 


# Citation
For citation and further information please refer to the preprint:

Fan, B. et al. RNA-Puzzles Round V: Blind predictions of 23 RNA structures. Preprint.

Thank you for your interest in our work. We hope you find the code and resources provided here useful in reproducing and building upon the findings of the RNA-Puzzles Round V.

Hajdin et al., RNA (7) 16, 2010  
RNA. 2009 Oct; 15(10): 1875–1885.

# Issues
If you encounter any issues or have questions about the code or analysis, please open an issue on the GitHub repository.
