# butyrateGutBoneAxis
Three compartment pharmacokinetic model of butyrate in the gut-bone axis

[![DOI](https://zenodo.org/badge/384794234.svg)](https://zenodo.org/badge/latestdoi/384794234)

## Overview 
This mathematical model describes the effect of butyrate on regulatory T cells in the gut, blood, and bone compartments. This model also explores the direct and indirect immune-mediated impacts of butyrate on bone metabolism via TGF-β and Wnt10b signaling molecules.

## Mathematical modeling of the gut-bone axis and implications of butyrate treatment on osteoimmunology
### Code Authors
Mohammad Aminul Islam and Ashlee N. Ford Versypt, 
Dept. of Chemical & Biological Engineering,
University at Buffalo, The State University of New York.

Corresponding author: A. N. Ford Versypt, ashleefv@buffalo.edu

### Manuscript
M. A. Islam, C. V. Cook, B. J. Smith, and A. N. Ford Versypt, Mathematical Modeling of the Gut-Bone Axis and Implications of Butyrate Treatment on Osteoimmunology, Industrial & Engineering Chemistry Research, 2021,60(49), 17814–17825. [DOI: 10.1021/acs.iecr.1c02949.](https://pubs.acs.org/doi/abs/10.1021/acs.iecr.1c02949)

### Scripts

* butyrate_model.py
This file contains the function for obtaining the concentration profile of butyrate, percentage change in regulatory T cells, fold change of TGF-β and Wnt-10b, and dynamics of bone cells and bone volume. Figures 3, 4, 5, 6, 11, and 12 are generated using this file. 

* butyrate_model_constants.py
This file contains input values and model parameters for butyrate_model.py file. The input value for change in intestine butyrate concentration can be varied from up to 0.2 μM.

* butyrate_fitting.py
This file contains the function for curve fitting with in vitro experimental results. Figure 2 is generated using this file.

* butyrate_model_curve_fitting_TGF.py
This file contains the function for simultaneous iterative curve fitting to estimate the dynamics of TGF-β for both positive and negative change of intestine butyrate concentration. Figure 10 is generated using this file.

* butyrate_optimizing_Wnt.py
This file contains the function for estimating the parameter ρ<sub>w </sub> for Wnt-10b.
  
* butyrate_optimizing_Wnt_constants.py
This file contain corresponding fixed parameters for butyrate_optimizing_Wnt.py.

* local_Sensitivity.py
This file contains the function for local sensitivity analysis. Figure 7 is generated using this file.

* global_Sensitivity.py
This file contains the function for global sensitivity analysis. Figure 8 is generated using this file.

* prediction_Tregs.py
This file contains the function for predicting change in bone regulatory T cells for any change in intestine butyrate concentration. Figure 9 is generated using this file.

* prediction_Tregs_constants.py
This file contains input value and model parameters for prediction_Tregs.py file.

## Acknowledgements
Research reported in this publication was supported by the National Institute of General Medical Sciences of the National Institutes of Health under award number R35GM133763. The content is solely the responsibility of the authors and does not necessarily represent the offcial views of the National Institutes of Health.
