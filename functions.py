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
        updated :
         n : 3D array (length, width, 9)
             shifted denisty array according to lattice vector 
             
    """
    
    for i in range(len(self.e)):
        
        par.n[:,:,i] = np.roll(par.n[:,:,i], self.e[i], axis = [0, 1])
        
    return par

def boundary_bounch(self, par):
    """ Takes the upper boundary, lower boundary and possible obstacle of a certain density vector. 
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
        updated :
        n : 3D array (length, width, 9)
            inverse directions for denisties in the boundary
    
    """  
    
    par.n[:,(0,-1),1:] = np.roll(par.n[:,(0,-1),1:], 4, axis = 2)
    
    if self.obs != 'none':
        par.n[par.indices,1:] = np.roll(par.n[par.indices,1:], 4, axis = 1)
    
    return par

def obstruction(self, par):
    """ Determines the weighted velocities for all directions for every grid point
    
    Parameters
    ----------
    self : NameSpace
        simulation constants
        obs : Str
            can be none, square or circle 
    par : NameSpace
        containing parameter arrays
        
    Returns
    -------
    par : NameSpace
        updated : 
        indices : 2D array
            contains the indices/coordinates of the obstacle

        
    """

    if self.obs == 'square':
        R = self.R
        cx, cy = int(self.L/4), int(self.W/2)

        X,Y = np.meshgrid(np.linspace(0, self.L, self.L_n), np.linspace(0, self.W, self.W_n))
        grid = np.stack((X,Y), axis = -1)
        grid[:,:,0] -= cx
        grid[:,:,1] -= cy

        par.indices = np.transpose(abs(grid[:,:,0]) + abs(grid[:,:,1]) <= R)
            
    elif self.obs == 'cylinder':
        R = self.R
        cx, cy = int(self.L/4), int(self.W/2)

        X,Y = np.meshgrid(np.linspace(0, self.L, self.L_n), np.linspace(0, self.W, self.W_n))
        grid = np.stack((X,Y), axis = -1)
        grid[:,:,0] -= cx
        grid[:,:,1] -= cy

        par.indices = np.transpose((grid[:,:,0]**2 + grid[:,:,1]**2) <= R**2)
    
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
        updated : 
        n_eq : 3D array (length, width, 9)
            contains new equilibrium distribution based on speed and density
        
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
        updated : 
        rho : 2D array (length, width)
            density at every grid point
        u : 3D array (length, width, 2)
            velocity for grid point in x [:,:,1] and y [:,:,0] direction
            weighted by density
        
    """

    par.rho = np.sum(par.n, axis = 2)
    par.u = np.tensordot(par.n, self.e, axes = 1)  
    
    valid_rho = par.rho != 0
    
    for k in range(len(par.u[0,0,:])):    
        par.u[valid_rho,k] = par.u[valid_rho,k]/par.rho[valid_rho]
    
    if self.obs != 'none':
        par.u[par.indices,:] = 0
        
    # Store velocity and density     
    par.v_tot.append(np.sum(abs(par.u))/(self.L*self.W))
    par.rho_tot.append(np.sum(par.rho)/(self.L*self.W))
    
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
        upated : 
        u : 3D array (length, width, 2)
            added forcing in the y direction [:,:,0]
        
    """
    
    par.u[:,1:self.W_in,0] = par.u[:,1:self.W_in,0] + self.dv*self.c
    
    if self.obs != 'none':
        par.u[par.indices, 0] -= self.dv*self.c
    
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
        updated :
        n : 3D array (lenth, width, 9)
            denisty relaxed to equilibrium denisty based on 
            relaxation time. 
        
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

