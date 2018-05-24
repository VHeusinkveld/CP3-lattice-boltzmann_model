import numpy as np

# -----------------------------------------------------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------------------------------------------------

def shift_n(self, par):
    """ Shifts densities according to their unit vectors. 
    
    Parameters
    ----------
    self : NameSpace
        simulation constants 
    par : NameSpace
        simulation parameters 
    
    Returns
    -------
    par : NameSpace
        containing parameter arrays and
        updated simulation density (n) parameters
    
    """
    
    for i in range(len(self.e)):
        
        par.n[:,:,i] = np.roll(par.n[:,:,i], self.e[i], axis = [0, 1])
        
    return par

def boundary_bounch(self, par):
    """
    
    Takes the upper boundary and the lower boundary of a certain density vector. 
    Mirrors this vector and assigns the density accordingly. 
    
    Parameters
    ----------
    self : NameSpace
        simulation constants 
    par : NameSpace
        simulation parameters 
    
    Returns
    -------
    par : NameSpace
        containing parameter arrays and
        updated simulation density (n) parameters
    
    """    
    
    for i, j in self.bounch:
        
        # Select upper and lower boundary 
        bd_1 = par.n[:,[0, self.W_in], i]
        bd_2 = par.n[:,[0, self.W_in], j]      
            
        # Exchange densities accordingly 
        par.n[:,[0, self.W_in], i] = bd_2
        par.n[:,[0, self.W_in], j] = bd_1
        
        # in progress for object 
        #bd_3 = par.n[par.obs, i]
        #bd_4 = par.n[par.obs, j] 
        #par.n[par.obs, i] = bd_4
        #par.n[par.obs, j] = bd_3
    
    
    return par

def obstruction(self, par):
    """ WIP function concerning an object in the pipe.
    
    """
    a = 4
    b = 4

    x = np.arange(30, 30 + int(self.L_in/a))
    y = np.arange(10, 10 + int(self.W_in/b))
    
    par.obs = np.meshgrid(x, y)
    
    return par

def eq_n(self, par):
    """ Determines the equilibrium densities for all directions in grid
    
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
        updated simulation equilibrium density (n_eq) and velocity (u) parameters
        
    """
    par.n_eq = np.zeros(np.shape(par.n), dtype = float)
    
    c = self.c
    u = par.u
 
    for i in range(len(self.e)):
        
        par.n_eq[:,:,i] = self.w[i]*par.rho/self.m*(1 + (3/c**2)*np.dot(u, self.e[i]) +
                                                    (9/(2*c**4))*np.dot(u, self.e[i])**2 -
                                                    (3/(2*c**2))*np.sum((u * u), axis = 2))                      
    return par 

def velocity(self, par):
    """ Determines the weighted velocities for all directions for every grid point
    
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
        updated simulation velocity (u) and average density (rho) parameters
        
    """

    par.norm_v = np.tensordot(par.n, self.e_norm, axes = 1)
    par.u = np.tensordot(par.n, self.e, axes = 1)  
    
    valid_rho = par.norm_v != 0
    #invalid_rho = par.norm_v == 0
    
    # Try to remove this for loop
    for k in range(len(par.u[0,0,:])):
        
        par.u[valid_rho,k] = par.u[valid_rho,k]/par.norm_v[valid_rho]
     
    par.v_tot.append(np.sum(abs(par.u))/(self.L*self.W))
    
    par.rho = np.sum(par.n, axis = 2) 
    
    return par

def forcing(self, par):    
    """ Changes the velocity according to the pressure gradient in the system
    
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
        updated simulation velocity (u) parameters
        
    """
    
    par.u[:,1:self.W_in,0] = par.u[:,1:self.W_in,0] + self.dv*self.c

    return par 

def relax_n(self, par):
    """ Determines the densities in the system after relaxation
    
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
        updated simulation density (n) parameters
        
    """
     
    par.n = (1 - 1/self.tau)*par.n + par.n_eq/self.tau
    
    return par