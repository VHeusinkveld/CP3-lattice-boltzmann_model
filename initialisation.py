import numpy as np
from sys import exit 
from types import SimpleNamespace
from functions import forcing 
from functions import eq_n
from functions import obstruction

# -----------------------------------------------------------------------------------------------------------------------
# Input parameter checks
# -----------------------------------------------------------------------------------------------------------------------

def input_check(self):
    """ Checks imput constants and gives an error when 
        - Width, height or res have wrong values. 
        - The relaxation time is shorter then a simulation time step.
        
    Parameters
    ----------
    self : NameSpace
        simulation constants/properties 
        
    Returns
    -------
    error warning if needed
        
    """
    
    if (self.L/self.res)%1 !=0 or (self.W/self.res)%1 != 0:
        exit('Choose width, height and res such that an integer ammount of points is generated.')
        
    if self.tau <= 1:
        exit('Tau should be bigger then 1.')

# -----------------------------------------------------------------------------------------------------------------------
# Initialisation functions
# -----------------------------------------------------------------------------------------------------------------------

def initialization(self):
    """ Initializes arrays to initial conditions. 
    Density n is set to the equilibrium density, n_eq.
    Borders and obstacles are initialised with 0 denisty .
    Velocity u starts at 0 but a forcing is applied once in the y direction.
    
    Parameters
    ----------
    self : NameSpace
        simulation constants 
        
    Returns
    -------
    self : NameSpace
        updated : 
        grid : 2D array (length, width)
            contains grid coordinates 
        obs : 2D array (length, width)
            contains obstacle coordinates
    
    par : NameSpace 
        
        updated :
        
        v_tot : list
            to keep track of speeds
        rho_tot : list
            to keep track of density
        kin_tot : list
            to keep track of kinetic energy   
        n : 3D array (length, width, 9)
            initialisation of the densities for the nine lattice vectors 
        u : 3D array (length, width, 2)
            x and y speed for every lattice point 
    
    
    """
    
    # Simulation parameters 
    par = SimpleNamespace()
    
    # Keeping track of conservation of quantities
    par.v_tot = []
    par.rho_tot = []
    par.kin_tot = []
    
    # Cartesian grid coordinates 
    self.grid_coord = np.meshgrid(np.linspace(0, self.L, self.L_n), np.linspace(0, self.W, self.W_n))
    
    # Integer grid coordinates
    self.grid_int_coord = np.meshgrid(range(self.L_n), range(self.W_n))
    
    # Density array
    par.n = np.zeros((self.L_n, self.W_n, len(self.e)), dtype = float)
    par.rho = np.ones ((self.L_n, self.W_n), dtype = float)
    par.u = np.zeros((self.L_n, self.W_n, 2), dtype = float)
    
    # Adds obstruction 
    if self.obs:
        par = obstruction(self, par) 
    
    # Initial conditions 
    par = forcing(self, par)    
    par = eq_n(self, par)
    par.n = par.n_eq
    
    if self.obs != 'none':
        par.rho[par.indices] = 0
        par.n[par.indices] = 0
    
    par.rho[:,[0, self.W_in]] = 0
    par.n[:,[0, self.W_in]] = 0
    
    return self, par