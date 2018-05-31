import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------------------------------------------------
# Plot functions
# -----------------------------------------------------------------------------------------------------------------------
  
def plot_markup():
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('font', size=18)

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
    plot_markup()
    
    plt.subplot(1, 3, 1)
    plt.imshow(par.rho)
    #plt.xlabel('ylabel')
    #plt.ylabel('xlabel')
    plt.title('rho')
    plt.colorbar()
    plt.subplot(1, 3, 2)
    plt.imshow(par.u[:, 1:self.W_in, 0])
    #plt.xlabel('ylabel')
    #plt.ylabel('xlabel')
    plt.title('vy')
    plt.colorbar()
    plt.subplot(1, 3, 3)
    plt.imshow(par.u[:, 1:self.W_in, 1])
    #plt.xlabel('ylabel')
    #plt.ylabel('xlabel')
    plt.title('vx')
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
    plot_markup()
    
    if self.obs == 'none':
        data_slice = par.u[0,1:self.W_in,0]
        plt.plot(data_slice, label = 'Simulation')

        x = np.linspace(0, len(data_slice) - 1, len(data_slice))
        y = -self.dv/(12*self.nu)*(x-(len(data_slice) - 1)/2)**2 
        y = y - min(y) + min(data_slice)
        plt.plot(x, y, '--', label = 'Theoretical')
        plt.xlabel('x')
        plt.ylabel('v_y')
        plt.legend()
        plt.tight_layout()
        if self.save:
            plt.savefig(self.fig_dir +'velocity_slice.png')
        plt.close()

    
    