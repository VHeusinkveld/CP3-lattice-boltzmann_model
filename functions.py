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
        
        par.n[:,:,i] = np.roll(par.n[:,:,i], self.e[i], axis = [1, 0])
        
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
    
    # To keep track of exchanged boundaries
    exchanged = np.ones((len(self.e),), dtype = bool)
    
    for i in range(len(self.e)):
        if exchanged[i]:
            e_inv = -1 * self.e[i]
            j = np.argwhere(np.sum(e_inv == self.e, axis = 1) == 2)[0]
            
            # Select upper and lower boundary 
            bd_1 = par.n[:,[0, self.W_in], i]
            bd_2 = par.n[:,[0, self.W_in], j]
            
            # Exchange densities accordingly 
            par.n[:,[0, self.W_in], i] = bd_2
            par.n[:,[0, self.W_in], j] = bd_1
            
            exchanged[i] = False
            exchanged[j] = False
            
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
      
    par.rho = np.tensordot(par.n, self.e_norm, axes = 1)
    
    par.u = np.tensordot(par.n, self.e, axes = 1)
    
    # Try to remove this for loop
    for k in range(len(par.u[0,0,:])):

        par.u[:,:,k] = par.u[:,:,k]/par.rho

    return par

def forcing(self, par):
    #incorrect still, not sure if this is a x or y forcing 
    
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
    
    par.u[:,:,1] = par.u[:,:,1] + 0.5*self.c #(Part of 0.5*self.c should be a relation with pressure instead. See ex. Ch. 14)

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