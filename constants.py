# -----------------------------------------------------------------------------------------------------------------------
# System constants
# -----------------------------------------------------------------------------------------------------------------------
import numpy as np
from types import SimpleNamespace

def constants(self):
    """ Sets the system constants based on the set system parameters
    
    Parameters
    ----------
    self : NameSpace
        set simulation constants
        
    Returns
    -------
    self : NameSpace
        updated simulation constants
        
    """
    # Resolutions
    self.res = 1 # Space resolution
    self.dt = 1 # Time resolution

    # Grid constants
    self.L_in, self.W_in = int(self.L/self.res), int(self.W/self.res)  
    self.L_n,  self.W_n  = self.L_in + 1, self.W_in + 1

    # Weights and direction vectors
    #self.w = np.array([    4,   1/4,     1,   1/4,      1,    1/4,       1,    1/4,      1])/9 # Jos vectors
    self.w = np.array([    4,     1,   1/4,     1,    1/4,      1,     1/4,      1,    1/4])/9  # Online source
    self.e = np.array([[0,0], [1,0], [1,1], [0,1], [-1,1], [-1,0], [-1,-1], [0,-1], [1,-1]])
    self.e_norm = np.sum(abs(self.e), axis = 1)
    self.e_norm = np.ones(np.shape(self.e_norm)) #to set the norm to 1
    self.bounch = np.array([[1,5], [2,6], [3,7], [4,8]])

    self.c = self.res/self.dt 
    self.tau = (6*(self.nu*(self.dt/self.res**2)) + 1)/2 # Relaxation time
    
    return self