# pl-validation
Verifying that persistence length calculations of conjugated donor acceptor polymers using flowermd are independent of each other.  
This repository is based on work completed by Dr. Paul found here: https://github.com/cmelab/forcefields/tree/main  
Polymers being studied (Danielsen et al.): https://pubs.acs.org/doi/10.1021/acs.macromol.1c02229#_i2  

The ff_generation.ipynb generates mbuild monomers and their forcefields using [ESPALOMA](https://docs.espaloma.org/en/latest/). A foyer forcefield is output that should work with an uncharged mol2 file. 

poly_builder.ipynb uses mbuild to polymerize the monomers.

poly_min_working_example.ipynb contains code to run NVT simulations.

coarse_grain.ipynb makes CG beads in preparation for persistence length calculations.

data_analysis.ipynb has the requisite functions for calculating the persistence length of mbuild polymers.

## To install environment:
```
git clone https://github.com/cmelab/pl-validation.git
conda env create -f environment.yml
conda activate p_l
```
