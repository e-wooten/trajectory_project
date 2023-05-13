This code computes the kinetic energy (kj/mol), potential energy (kj/mol), total energy (kj/mol), temperature (K), pressure (atm), volume (A^3), and root mean square deviation (RMSD) angstroms over the course of a molecular dynamic simulation. Plots are made to show each of the properties mentioned above over the course of a molecular dynamic simulation.

Author: Eric Wooten
Date: 04/2022

There are 3 input files needed to run this program.

1) ".log" file
A log file lists the parameters of a molecular simulation. It takes data such as pdb files and psi files and uses the parameters for the simulation. It outputs specified time steps and other pertinent information such as energies, pressure, volume, temperature, bonds, angles, dihedrals. Data is written to this file in time steps with each of the previously mentioned properties at each time step. 

2) ".psf" file
A psf is a protein structure file. A protein structure files contains topology information for force fields of the protein structure. Included in this is bonding information that is used in molecular simulations in order for certain parameters to be defined. Some information in a psf is connectivity and charge.

3) ".dcd" file
A dcd file is also known as a trajectory file. It is basically a pdb for every snapshot throughout the molecular simulation. It is one of the many types of trajectory files. This file is what allows us to visualize the simulation.

Usage:
To use this code you simply need to supply a path to a .log file from a simulation as well as paths to .psf and .dcd files. Once these files are input into the program the out put will display 7 graphs with kinetic energy (kj/mol), potential energy (kj/mol), total energy (kj/mol), temperature (K), pressure (atm), volume (A^3), and root mean square deviation (RMSD) angstroms all plotted against the time (step) of the simulation.


There are 5 necessary libraries that must be imported at the start to execute the program.
 
1)Import os
-This library allows us to define the path to the file in the directory it is in.

2)import pandas as pd
- Pandas is a data manipulation tool. It easily stores data making it accessible to all sorts of calculations and plotting

3) import MDAnalysis as mda
_This library is useful for storing trajectory files along with a psf and get meaningful information from the trajectory. It can retrieve things such as information of the atoms of the system, residues in the protein, as well properties such as charge. 

4)from MDAnalysis.analysis import rms
-This library allows the importation of a root mean square deviation (RMSD) calculator over the course of a trajectory (dcd) file

5) import matplotlib.pyplot as plt
-This library is essential for plotting all of the information received from the input files.

Initial input:
The initial inputs are a .log, .psf, and .dcd file. There are 3 examples below of the files of lysozyme

Executing the code:
-To execute the code all that is needed to be done is enter the specific .log, .psf., and .dcd files as prompted.

- Once the log file is input it is read and becomes a (data) variable.
- The (data) variable is then iterated through to find time (step), kinetic energy (kj/mol), potential energy (kj/mol), total energy (kj/mol), temperature (K), pressure (atm), and volume (A^3). These outputs all become their own variable lists in each list is the attributed value of each time step im the log file.
-These variable lists for each property are then zipped into a pandas data frame. This outputs a data frame of all these property values over the entire course of the simulation.
-Once in the data frame, the information is used to output plots of kinetic energy (kj/mol), potential energy (kj/mol), total energy (kj/mol), temperature (K), pressure (atm), and volume (A^3) vs time (step) of the simulation.   

-The RMSD output is gathered from the .psf and .dcd files. Once the .psf and .dcd files are input the files are read by the MDAnalysis library and saved as the variable (protein_simulation)
- This is then input into the rmsd calculation. Each structure of the trajectory is aligned by alpha carbons with the reference frame, in this case frame 0. Other selections can also have the RMSD calculated such as 'protein backbone' and 'all atoms of the protein'.
- This code selects alpha carbons, protein backbone, and all protein atoms to make an RMSD calculation over the course of the trajectory.
- These RMSD values over the course simulation are then stored in a pandas data frame consisting of columns [Frame number, Time (ps), Alpha Carbons, Protein backbone atoms, and all protein atoms.
- Matplot lib is then used to plot the RMSD of the three choices above in the same group. 


Example:

Example Inputs: (These inputs are provided to run the program)
Enter the path to the log file that contains the data: /Users/kingstreasure/Documents/programming_for_data_analysis/trajectory_project/lza_wb_eq.log

Enter the path to the psf file that contains the data: /Users/kingstreasure/Documents/programming_for_data_analysis/trajectory_project/lza_wb_solvate.psf

Enter the path to the dcd file that contains the data: /Users/kingstreasure/Documents/programming_for_data_analysis/trajectory_project/lza_wb_eq.dcd

Outputs: (The outputs are then plots of the information found in these example output files)
log_output.png
rmsd_output.png
    




# trajectory_project
# trajectory_project
