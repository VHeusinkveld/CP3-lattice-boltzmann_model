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
    par.n[:,(0,-1),1:] = np.roll(par.n[:,(0,-1),1:], 4, axis = 2)
    '''
    for i, j in self.bounch:
        
        # Select upper and lower boundary 
        bd_1 = par.n[:,[0, self.W_in], i]
        bd_2 = par.n[:,[0, self.W_in], j]      
            
        # Exchange densities accordingly 
        par.n[:,[0, self.W_in], i] = bd_2
        par.n[:,[0, self.W_in], j] = bd_1
        
        # in progress for object 
        if self.obs:
            bd_3 = par.n[par.obs_x, par.obs_y, i]
            bd_4 = par.n[par.obs_x, par.obs_y, j] 
            par.n[par.obs_x, par.obs_y, i] = bd_4
            par.n[par.obs_x, par.obs_y, j] = bd_3
    '''
    
    return par

def obstruction(self, par):
    """ WIP function concerning an object in the pipe.
    
    """
   
    x_start = int(self.L_in/4 - 1/2*self.L_in*self.L_ratio)
    y_start = int(self.W_in/2 - 1/2*self.W_in*self.W_ratio)
    
    x = np.arange(x_start, x_start + int(self.L_in*self.L_ratio))
    y = np.arange(y_start, y_start + int(self.W_in*self.W_ratio))
    
    par.obs = np.meshgrid(x, y)
    par.obs_x = np.reshape(par.obs[0],(-1,))
    par.obs_y = np.reshape(par.obs[1],(-1,))
    
    x_int = np.arange(x_start + 1, x_start - 1 + int(self.L_in*self.L_ratio))
    y_int = np.arange(y_start + 1, y_start - 1 + int(self.W_in*self.W_ratio))
    
    par.obs_int = np.meshgrid(x_int, y_int)
    par.obs_int_x = np.reshape(par.obs_int[0],(-1,))
    par.obs_int_y = np.reshape(par.obs_int[1],(-1,))
    
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
    
    for k in range(len(par.u[0,0,:])):
        
        par.u[valid_rho,k] = par.u[valid_rho,k]/par.norm_v[valid_rho]
    
    if self.obs:
        par.u[par.obs_int_x, par.obs_int_y,:] = 0
        
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
    
    if self.obs:
        par.u[par.obs_x, par.obs_y, 0] += -self.dv*self.c
    
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

def Re_pipe(self, par):   
    return np.mean(par.u[:,:,0])*self.W/self.nu

def Reynolds(self, par):
    Reynolds_pipe = Re_pipe(self, par)
    print('----------------------------')
    print('Reynolds_pipe = ' +str(Reynolds_pipe))
    
    if self.obs:
        Reynolds_obs = Re_obs(self, par)
        print('Reynolds_obs = ' +str(Reynolds_obs))
        
    else:
        Reynolds_obs = np.nan
        
    return Reynolds_pipe, Reynolds_obs

# -----------------------------------------------------------------------------------------------------------------------
# Under development
# -----------------------------------------------------------------------------------------------------------------------

def Re_obs(self, par):
    return np.mean(par.u[0:(par.obs_x[0]-5), par.obs_y, 0])*(self.W_ratio*self.W)/self.nu

