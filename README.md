# MassSpringModels_ovuleGrowth2D

[![DOI](https://zenodo.org/badge/233707711.svg)](https://zenodo.org/badge/latestdoi/233707711)

plugin for MorphoMechanX(www.MorphoMechanX.org) to generate 2D simulations of growing and dividing ovules at cellular resolution based on mass-springs.

To run this model on the test mesh provided, get MorphoMechanX first (you can find it here: www.morphomechanx.org). As a reference, this model is guaranteed to run against version 2.0.2-858 (git commit: 8eb2aec9055221e3e351ffe064bd059186717374)

The Mass-Spring model itself does not require installation.
After MorphoMechanX has been installed, to run the model:

1) open a terminal from the "model" directory provided in this repository

2) type "make clean" from the terminal

3) type "make run" from the terminal.

This will open a graphically interactive window of MorphoDynamX, and MorphoMechanX tab will be shown (under Process/Model).

Depending on which specific model the user wants to run (with reference to the model variants defined in the paper), choose a model number (for example "MS-Model2"):

-open the chosen model folder from a navigation browser (i.e. a terminal from "MS-Model2")

-drag and drop the view file: "CellDisk.mdxv"  into the already open MorphoDynamX window 

(altrernatively go into MorphoDynamX and open the File tab and from there select "CellDisk.mdxv" from the chosen folder case).

By double-clicking on: 

00 Ovule Master Process (Process/Model/ Cell Ovule Growth/00 Ouvle Master Process)
the simulation will start selecting the right template from the model variant chosen. 

For further details on how to edit model parameters or their meaning, refer to the Appendix 1 provided with the paper (Hernandez-Lagana, Mosca, Mendocilla-Sato et al.). The Appendix can be found as pre-print also here https://www.biorxiv.org/content/10.1101/2020.07.30.226670v3.supplementary-material under "Supplemental File 1"

For further questions, contact: gabriella.mosca@gmail.com
