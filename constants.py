# -----------------------------------------------------------------------------------------------------------------------
# System constants
# -----------------------------------------------------------------------------------------------------------------------
from types import SimpleNamespace

# Resolutions
sim.res = 1 # Space resolution
sim.dt = 1 # Time resolution

# Grid constants
sim.L_in, sim.W_in = int(sim.L/sim.res), int(sim.W/sim.res)  
sim.L_n,  sim.W_n  = sim.L_in + 1, sim.W_in + 1

# Weights and direction vectors
sim.w = np.array([4, 1/4, 1, 1/4, 1, 1/4, 1, 1/4, 1])/9
sim.e = np.array([[0,0], [1,0], [1,1], [0,1], [-1,1], [-1,0], [-1,-1], [0,-1], [1,-1]])
sim.e_norm = np.sum(abs(sim.e), axis = 1) # takes care of normalization since some vectors are [1, 1]

sim.c = sim.res/sim.dt 
sim.tau = (6*(sim.nu*(sim.dt/sim.res**2)) + 1)/2 # Relaxation time