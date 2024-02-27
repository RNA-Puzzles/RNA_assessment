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
This package depends on `biopython`.     
To install: `pip install biopython`   
OpenStructure for lDDT calculations.

If need to calculate the Interaction Network Fidelity, it needs to call [`MC-annotate`](https://major.iric.ca/MajorLabEn/MC-Tools.html).    
Please download the binary execution from the website and coordinate the directory for it at the top line `MCAnnotate_bin=` of the mcannotate.py script. 

If need to calculate the Local Distance Difference Test (lDDT), it needs to call OpenStructure (https://openstructure.org/download).    

If need to calculate the Atomic Rotationally Equivariant Scorer (ARES), it needs to call ARES (https://www.science.org/doi/10.1126/science.abe5650)
        
        
        


## how to use
A detailed introduction can be found in the example [notebook](https://github.com/RNA-Puzzles/RNA_assessment/blob/master/example.ipynb) or the example [script](https://github.com/RNA-Puzzles/RNA_assessment/blob/master/example/example.py). 

## citation
Hajdin et al., RNA (7) 16, 2010  
RNA. 2009 Oct; 15(10): 1875–1885.
