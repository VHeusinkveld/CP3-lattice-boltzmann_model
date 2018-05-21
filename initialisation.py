import numpy as np
from sys import exit 
from types import SimpleNamespace

# -----------------------------------------------------------------------------------------------------------------------
# Input parameter checks
# -----------------------------------------------------------------------------------------------------------------------

def input_check(self):
    """ Checks imput constants and gives an error when width, height or res have wrong values
    
    Parameters
    ----------
    self : NameSpace
        simulation constants 
        
    Returns
    -------
    error warning if needed
        
    """
    
    if (self.L/self.res)%1 !=0 or (self.W/self.res)%1 != 0:
        exit('Choose width, height and res such that an integer ammount of points is generated.')

# -----------------------------------------------------------------------------------------------------------------------
# Initialisation functions
# -----------------------------------------------------------------------------------------------------------------------

def initialization(self):
    """ Initializes arrays. 
    
    Parameters
    ----------
    self : NameSpace
        simulation constants 
        
    Returns
    -------
    par : NameSpace
        containing initialized parameter arrays
    
    """
    
    # Simulation parameters 
    par = SimpleNamespace()
    
    # Cartesian grid coordinates 
    self.grid_coord = np.meshgrid(np.linspace(0, self.L, self.L_n), np.linspace(0, self.W, self.W_n))
    
    # Integer grid coordinates
    self.grid_int_coord = np.meshgrid(range(self.L_n), range(self.W_n))
    
    # Density array
    par.n = np.ones((self.L_n, self.W_n, len(self.e)), dtype = float)
    par.n[:,[0, self.W_in],:] = 0
    
    return self, par