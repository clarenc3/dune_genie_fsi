# dune_genie_fsi
GENIE FSI Reweighting scheme for DUNE TDR

Shifts the Eav/q0 spectrum for one q0, q3 bin by looking at the effects of FSI on and FSI off cross-sections.

The scheme maintains the total cross-section and shifts the mean and variance of the spectrum to a user-specified number (`meanshift` and `varshift`).

Uses a CC-inclusive selection, excluding CC-coherent.

Depends on [ROOT](http://root.cern.ch) and takes [NUISANCE](http://nuisance.hepforge.org) flat-trees as input.

## fsianal.cpp
fsianal.cpp produces differential cross-sections in q3, q0, Eav/q0 and auxillary histograms needed for the method. A "helpful" script (runit.sh) is provided to run on multiple files.
fsianal.cpp is run by e.g. `$ root -l -b -q 'fsianal.cpp("NUISANCE_FLAT_TREE.root")'`

## KevinRecipe.cpp
Then KevinRecipe.cpp takes the outputs of fsianal.cpp and produces the weighting function in q3, q0, Eav/q0 and some validation plots, including a multi-page pdf of the cross-section distributions in Eav/q0 for a given q3, q0 bin. 
KevinRecipe.cpp is run by e.g. `$ root -l -b -q 'KevinRecipe.cpp("12C_FSI_ON.root", "12C_FSI_OFF.root", "40Ar_FSI_ON.root")'`

[Talk to me](mailto:cwret@fnal.gov)
