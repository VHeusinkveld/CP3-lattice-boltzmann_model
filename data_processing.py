import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------------------------------------------------
# Plot functions
# -----------------------------------------------------------------------------------------------------------------------

def plot_boltzmann_lattice(self, par):
    """ Plot of 2D boltzman lattice with velocity profile over entire lattice
    
    Parameters
    ----------
    self : NameSpace
        simulation constants
    par : NameSpace
        containing parameter arrays
        
    Returns
    -------
    plots of system (UPDATE THIS DESCRIPTION!!!!!!!)
        
    """
    data_x = par.u[:, 1:self.W_in, 0]
    data_x[par.obs_x, par.obs_y] = 0

    plt.subplot(1, 2, 1)
    plt.imshow(data_x)
    plt.xlabel('ylabel')
    plt.ylabel('xlabel')
    plt.colorbar()
    plt.subplot(1, 2, 2)
    plt.imshow(par.u[:, 1:self.W_in, 1])
    plt.xlabel('ylabel')
    plt.ylabel('xlabel')
    plt.colorbar()
    plt.show()
    
def plot_velocity_profile(self, par):
    """ Plot of velocity profile over the width of the lattice
    
    Parameters
    ----------
    self : NameSpace
        simulation constants
    par : NameSpace
        containing parameter arrays
        
    Returns
    -------
    plots of system (UPDATE THIS DESCRIPTION!!!!!!!)
        
    """
    #x_axis = np.arange(self.W_in)
    plt.plot(par.u[0,1:self.W_in,0])
    plt.xlabel('Width')
    plt.ylabel('Velocity')
    plt.title('Velocity profile at t = ...')
    plt.show()