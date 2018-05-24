import numpy as np
from sys import exit 
from types import SimpleNamespace
import matplotlib.pyplot as plt

#from constants import *
from data_processing import *
from functions import *
from initialisation import *

# -----------------------------------------------------------------------------------------------------------------------
# Simulation
# -----------------------------------------------------------------------------------------------------------------------
    
def boltzmann_sim(self):
    """ Simulation of the boltzmann lattice
    
    Parameters
    ----------
    self : NameSpace
        simulation constants 
        
    Returns
    -------
    par : NameSpace
        containing parameter arrays in equilibrium stage
        
    """
    
    input_check(self)
    sim, par = initialization(self)
    
    counter = 0
        
    while True:
        counter += 1
        simulation_step(self, par)
        if counter == 1:
            par.u_reference = par.u*0 # used for equilibrium determination (Need array of same dimensions)
            plot_boltzmann_lattice(self, par)
            #plot_velocity_profile(self,par)
        par, equilibrium = is_stable(self, par)
        
        if counter%self.plot_iteration == 0:
            plot_boltzmann_lattice(self, par)
        if equilibrium:
            plot_boltzmann_lattice(self,par)
            #plot_velocity_profile(self,par)
            print('Equilibrium has been reached after ' + str(counter) + ' iterations.')
            return par
            #break #only needed when nothing there is no return command
        if counter == self.max_iterations:
            plot_boltzmann_lattice(self,par)
            #plot_velocity_profile(self,par)
            print('Maximum number of iterations (' + str(counter) + ') have been reached. It is adviced to increase the maximum number of iterations or increase the error tolerance (epsilon).')
            return par
            #break #only needed when nothing there is no return command

def is_stable(self, par):
    '''Determines whether system is in equilibrium or not by calculating the differences
    in velocity for every grid point between each iteration. Once this difference is smaller
    than epsiolon for every grid point, 'equilibrium' is set to True.
    
    Parameters
    ----------
    self : NameSpace
        simulation constants
    par : NameSpace
        containing parameter arrays
        
    Returns
    -------
    par : NameSpace
        containing parameter arrays and
        updated reference velocity for the next iteration (u_reference)
    equilibrium : Boolean
        set to True once system is in equilibrium. False otherwise
        
    '''
    du = par.u - par.u_reference
    equilibrium = np.all(du < self.epsilon)
    par.u_reference = par.u
    
    return par, equilibrium

def simulation_step(self, par):
    '''Runs one simulation iteration
    
    Parameters
    ----------
    self : NameSpace
        simulation constants
    par : NameSpace
        containing parameter arrays
        
    Returns
    -------
    par : NameSpace
        containing updated parameter arrays
        
    '''
    
    par = shift_n(self, par)
    par = boundary_bounch(self, par)
    par = velocity(self, par)
    par = forcing(self, par)
    par = eq_n(self, par)
    par = relax_n(self, par)
    
    return par




    