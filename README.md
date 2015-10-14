# UQOpenFOAM
Turbulence models for parametric UQ studies in OpenFOAM using the Bayesian-Model-Scenario-Averaging (BMSA) framework by
Edeling et al. (2014)

Features
-----------
The following turbulence models are available in the library

* Wilcox (2006) k-omega (Based on the implementation by Gomez et al., see: http://turbmodels.larc.nasa.gov/ChangesToOpenFOAM.pdf)
* Launder-Sharma Low-Re k-epsilon (Based on the standard implementation in OpenFOAM 2.3.3)
* Spalart-Allmaras (Based on the standard implementation in OpenFOAM 2.3.3)


Installation
-----------
Each model can be compiled with the standard procedure

```bash
source <path-to-your-OpenFOAM-installation>/etc/bashrc
cd <path-to-model-in-reposotory>
wmake
```




References
-----------
W. Edeling

S. Gomez
