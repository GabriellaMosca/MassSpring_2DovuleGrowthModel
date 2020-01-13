# MassSpringModels_ovuleGrowth2D
plugin for MorphoMechanX(http://www.mpipz.mpg.de/MorphoGraphX/MorphoMechanX/) to generate 2D simulations of growing and dividing ovules at cellular resolution based on masssprings.

The Git revision number these models are guarantedd to run with is: 
e618301477af8a0bda2c344a7312195ba5acf379

The Git revision number should be provided together with the request for MorphoMechanX, which is to be addressed to Richard.Smith@jic.ac.uk
Installation instructions can be requested together with the package or asked alternatively to this email address: gabriella.mosca@gmail.com 

MorphoMechanX work on Linux and the models have been tested on Linux Mint 19. 

The plugin itself does not require installation. After MorphoMechanX has been installed, to run the model 
1)move into the model directory from a terminal
2)type "make run" and enter.  
This will open a graphically interactive window of MorphoDynamX, and MorphoMechanX tab will be shown (Process/Model).
Depending on which model case the user wants to run, choose a case number (i.e. case2): 
-open a windown into its folder (i.e. a window inside case2),
-drag and drop the view file: "CellDisk.mdxv" into MorphoDynamX window 
(altrernatively go into MorphoDynamX and open the File tab and from there select "CellDisk.mdxv" from the chosen folder case.
By double-clicking on: 
00 Ovule Master Process (Cell Ovule Growth/00 Ouvle Master Process)
the simulation will start selecting the right template from the case chosen. 

