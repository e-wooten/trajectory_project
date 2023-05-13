#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 03:20:40 2023

@author: kingstreasure
"""

#Importing necessary libraries to run the code
import os
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt


# Using the abspath function to define file location using absolute path and reading the log file
# Setting all the contents of the log file to the variable data

def read_log_file(path_to_log_file):
    """
    Reads the path_to_file of the log file and outputs the content in a data list

    INPUT
    ----------
    path_to_log_file: string
        Path to log file containing the data.

    OUPUT
    -------
    data: list
    list with the data in the log file

    """
    try:
        simulation_file_location = os.path.abspath(path_to_log_file)
        with open(simulation_file_location,'r') as simulation_file:
            data = simulation_file.readlines()
            simulation_file.close()
            return(data)
    except:
        print('Error reading path to file')

# Setting up empty lists of the relevant information we'd like to extract from the log file        
# Using a for loop to go through the log file and extract all information from that data that begins with the line ENERGY
# Splitting the line where energy was found so that every piece on information is its own string
# Changing the selected data from the split energy line into a float
# Appending the newly acquired information to the appropriate list
# Using the zip function to merge the lists together
# Transforming our zipped list into a dataframe with appropriate colum headings and units

def data_extraction(data):
    """
    Extracting timestep, potentital energy, kinetic energy, total energy, temperature, pressure, and volume from data into dataframe

    INPUT
    ----------
    data : list
        data is in list format.

    OUTPUT
    -------
    simulation_df: pandas dataframe
        simulation_df in pandas dataframe

    """
    try:
        time_step_list = []
        kinetic_energy_list = []
        potential_energy_list = []
        total_energy_list = []
        temperature_list = []
        pressure_list = []
        volume_list = []
        for line in data:
            if 'ENERGY:' in line:
                energy_line = line
                energy_line_split = energy_line.split()
                time_step = float(energy_line_split[1])
                time_step_list.append(time_step)
                kinetic_energy = float(energy_line_split[10])*4.184
                kinetic_energy_list.append(kinetic_energy)
                potential_energy = float(energy_line_split[13])*4.184
                potential_energy_list.append(potential_energy)
                total_energy = float(energy_line_split[11])*4.184
                total_energy_list.append(total_energy)
                temperature = float(energy_line_split[12])
                temperature_list.append(temperature)
                pressure = float(energy_line_split[16])/1.01325
                pressure_list.append(pressure)
                volume = float(energy_line_split[18])
                volume_list.append(volume)
        simulation_data = list(zip(time_step_list, kinetic_energy_list, potential_energy_list, total_energy_list, temperature_list,
                           pressure_list, volume_list))
        simulation_df = pd.DataFrame(simulation_data, columns = ['Time (step)', 'Kinetic Energy (kj/mol)', 'Potential Energy (kj/mol)',
                                                        'Total Energy (kj/mol)', 'Temperature (K)', 'Pressure (atm)', 'Volume (A^3)'])
        return(simulation_df)
    except:
        print('Error producing the simulation dataframe')
        
# Getting x,y1,y2,y3,y4,y5,y6 data from the simulation dataframe
# Making the figure with 2 subplots and specifying size and distance between plots
# Specifying the line plot, which x and y data that will be in each subplot, along with color and chart labels 

def plot_data(simulation_df):
    """
    Plotting time Kinetic Energy, Potential Energy, Total Energy (kj/mol), Temperature (K), Pressure (atm), Volume (A^3) vs Time (step) in individual plots from the simulation_df

    INPUT
    ----------
    simulation_df : pandas dataframe
        simulation_df in pandas dataframe.

    OUTPUT
    -------
    fig
        Subplots of kinetic energy, potential energy, total energy, temperature, pressure and volume over the entire simulation

    """
    try:
        x = simulation_df['Time (step)']
        y1 = simulation_df['Kinetic Energy (kj/mol)']
        y2 = simulation_df['Potential Energy (kj/mol)']
        y3 = simulation_df['Total Energy (kj/mol)']
        y4 = simulation_df['Temperature (K)']
        y5 = simulation_df['Pressure (atm)']
        y6 = simulation_df['Volume (A^3)']

        fig, ((ax1, ax2), (ax3, ax4),(ax5,ax6)) = plt.subplots(ncols=2, nrows=3,figsize=(10, 8))
        plt.subplots_adjust(wspace=0.3, hspace=0.4)

        ax1.plot(x, y1, '-', color='red')
        ax1.set_xlabel('Time (step)')
        ax1.set_ylabel('Kinetic Energy (kj/mol)')
        ax1.set_title('Kinetic Energy of the simulation')

        ax2.plot(x, y2, '-', color='blue')
        ax2.set_xlabel('Time (step)')
        ax2.set_ylabel('Potential Energy (kj/mol)')
        ax2.set_title('Potential Energy of the simulation')

        ax3.plot(x, y3, '-', color='purple')
        ax3.set_xlabel('Time (step)')
        ax3.set_ylabel('Total Energy (kj/mol)')
        ax3.set_title('Total Energy of the simulation')

        ax4.plot(x, y4, '-', color='gray')
        ax4.set_xlabel('Time (step)')
        ax4.set_ylabel('Temperature (K)')
        ax4.set_title('Temperature of the simulation')

        ax5.plot(x, y5, '-', color='green')
        ax5.set_xlabel('Time (step)')
        ax5.set_ylabel('Pressure (atm)')
        ax5.set_title('Pressure of the simulation')

        ax6.plot(x, y6, '-', color='orange')
        ax6.set_xlabel('Time (step)')
        ax6.set_ylabel('Volume ($\AA^3$)')
        ax6.set_title('Volume of the simulation')

        fig.savefig('/Users/kingstreasure/Documents/programming_for_data_analysis/project/log_output.png')
        fig.show()
    except:
        print('Error plotting the data')
        
# Using the abspath function to define file location using absolute path
# Saving the protein simulation structure and trajectories of the simulation to the variable 'protein_simulation'        
        
def read_simulation_files(path_to_psf_simulation_files, path_to_dcd_simulation_files):
    """
    Reads the path to the psf and dcd simulation files and outputs the content in a mda universe
    
    INPUT
    ----------
    path_to_simulation_files: string
        Path to the simulation files containing the data
        
    OUTPUT
    ----------
    protein_simulation: MDAnalysis.core.universe.Universe
        protein_simulation with number of atoms in the simulation is loaded
        
    """
    try:
        protein_simulation = mda.Universe(path_to_psf_simulation_files, path_to_dcd_simulation_files)
        return(protein_simulation)
    except:
        print('Error loading simulation files')
        
#Performs a rotational and translational alignment of the target trajectory to the reference structure
#Converting the results of the RMSD calculations between each frame and the reference frame into a pandas dataframe.
        
def rmsd_calculation(protein_simulation):
    """
    Calculates rmsd of the simulation for alpha carbons, protein backbone, and all protein atoms compared to the first frame.
    This is then put into a dataframe
    
    INPUT
    ----------
    protein_simulation: MDAnalysis.core.universe.Universe
        The protein_simulation data is stored as a universe in the MDA module
        
    OUTPUT
    ----------
    RMSD_df: pandas dataframe
        A data frame of calculated RMSD for alpha carbons, protein backbone, and all protein atoms compared to the first frame over the entire trajectory
        
    """
    try:
        RMSD = rms.RMSD(protein_simulation,  # structure that is being aligned
             protein_simulation,  # reference structure used for alignment
             #'Select' determines atoms are being used for alignemt, in this case (alpha carbons). Then RMSD is calculated for alpha carbons
             #'groupselections' are other criteria that RMSD is calculated, in this case backbone atoms and all protein atoms.
             select='name CA', groupselections=['backbone', 'protein'],  
             ref_frame=0)  # frame index of the reference
        RMSD.run()
        RMSD_df = pd.DataFrame(RMSD.results.rmsd,
                  columns=['Frame', 'Time (ps)', 'Alpha Carbons', 'Backbone atoms', 'All protein atoms'])
        return(RMSD_df)
    except:
        print('Could not calculate RMSD')
        
# Plotting the RMSD of the selection criteria against the frame number
# Specifying a line plot and labeling the plot
        
def rmsd_plot(RMSD_df):
    """
    Plots the calculated RMSD of the alpha carbons, protein backbone, and all protein atoms over the course of the simulation
    
    INPUT
    ----------
    RMSD_df: pandas dataframe
        A data frame of calculated RMSD for alpha carbons, protein backbone, and all protein atoms compared to the first frame over the entire trajectory
        
    OUTPUT
    ----------
    ax: plot
        A plot with the RMSD calculated values for alpha carbons, protein backbone, and all protein atoms over the course of the simulation
        
    """
    try:
        ax = RMSD_df.plot(x='Frame', y=['Alpha Carbons', 'Backbone atoms', 'All protein atoms'],
             kind='line')
        ax.set_ylabel(r'RMSD ($\AA$)');
        ax.set_title('RMSD throughout Simulation')
        plt.savefig('/Users/kingstreasure/Documents/programming_for_data_analysis/project/rmsd_output.png')
        plt.show();
        return(ax)
    except:
        print('Error plotting the RMSD data')
        

        
if __name__ == '__main__':
    try:
        
#Input data from the user
        path_to_log_file = input('Enter the path to the log file that contains the data')
               
#Read the log file
        data = read_log_file(path_to_log_file)
    
#Extract data from the log file
        simulation_df = data_extraction(data)
    
#Plot the data
        fig = plot_data(simulation_df)
    
#Input data from the user
        path_to_psf_simulation_files = input('Enter the path to the psf file')
        path_to_dcd_simulation_files = input('Enter the path to the dcd file')
        

#Read the simulation files
        protein_simulation = read_simulation_files(path_to_psf_simulation_files, path_to_dcd_simulation_files)
    
#Calculate RMSD values of the simulation
        RMSD_df = rmsd_calculation(protein_simulation)
               
#Plot the RMSD values
        ax = rmsd_plot(RMSD_df)


    except:
        print('Error running the program')