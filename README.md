# MassSpringModels_ovuleGrowth2D
plugin for MorphoMechanX(www.MorphoMechanX.org) to generate 2D simulations of growing and dividing ovules at cellular resolution based on mass-springs.

MorphoMechanX can be installed from package for Linux Mint 19 and Ubuntu 20.04. The package is provided at this git https://github.com/GabriellaMosca/PollenGrain_indentation (or https://zenodo.org/record/4590379#.YF3MyIBKiM8). The user should follow instructions from 1. to 4.

Alternatively, MMX can be installed from source code and can be requested at: www.MorphoMechanX.org.
The Git revision number these models are guarantedd to run with is: 

e618301477af8a0bda2c344a7312195ba5acf379

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

For detailed instructions on how to set the model parameters and their meaning, we refer to the "Appendix 1" provided with "Hernandez-Mosca-Mendocilla et al." 
 
