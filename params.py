import sys
import os
import logging
logging.root.setLevel(logging.INFO)
from distutils import spawn

import numpy as np
from crack_matscipy_update import G_to_strain
from ase.lattice import bulk
import ase.units as units
import ase.build
from ase.io import read

from qmmm import RescaledCalculator

# ********** Bulk unit cell ************

# 8-atom diamond cubic unit cell for diamond
a0_cast = 3.5694 # guess at lattice constant for carbon - we will minimize
a0_tersoff = 3.5
#cryst = bulk('C', 'diamond', a=a0, cubic=True)
cryst = ase.build.bulk('C', 'diamond', a=a0_cast, cubic=True)
surf_ny = 4 # number of cells needed to get accurate surface energy

#********** Diamond intrinisic parameters **********
E = 1281.1*units.GPa
nu = -0.01 


# *********  System parameters **********

# There are three possible crack systems, choose one and uncomment it

# System 1. (111)[0-11]
crack_direction = (-2, 1, 1)      # Miller index of x-axis
cleavage_plane = (1, 1, 1)        # Miller index of y-axis
crack_front = (0, 1, -1)          # Miller index of z-axis

# # System 2. (110)[001]
# crack_direction = (1,-1,0)
# cleavage_plane = (1,1,0)
# crack_front = (0,0,1)

# # System 3. (110)[1-10]
# crack_direction = (0,0,-1)
# cleavage_plane = (1,1,0)
# crack_front = (1,-1,0)

check_rotated_elastic_constants = False

width = 300.0*units.Ang              # Width of crack slab
height = 100.0*units.Ang              # Height of crack slab
depth = 2                            # The number of layers in the slab
vacuum = 50.0*units.Ang             # Amount of vacuum around slab
crack_seed_length = 150.0*units.Ang   # Length of seed crack
strain_ramp_length = 10.0*units.Ang  # Distance over which strain is ramped up
initial_G = 6.0*(units.J/units.m**2) # Initial energy flow to crack tip

relax_bulk = True                     # If True, relax initial bulk cell
bulk_fmax  = 1e-6*units.eV/units.Ang  # Max force for bulk, C_ij and surface energy

relax_slab = True                     # If True, relax notched slab with calculator
relax_fmax = 0.01*units.eV/units.Ang # Maximum force criteria for relaxation

# ******* Molecular dynamics parameters ***********

sim_T = 600.0*units.kB           # Simulation temperature
reset_temp = True
nsteps = 120                       # Total number of timesteps to run for
max_step_num = 1000               # Max number of iteration
timestep = 1.0*units.fs          # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang      # Amount by which potential cutoff is increased
                                 # for neighbour calculations
tip_move_tol = 10.0              # Distance tip has to move before crack
                                 # is taken to be running
strain_rate = G_to_strain(0.1, E, nu, height)


# ********* Defect parameters ************ 
defect = False 
vacancy = False 
defect_number = 7
defect_location = 14418

#*********** Writing parameters ***************
traj_interval = 1

# ******** Multiscale parameters *********

buffer_width = 7.0*units.Ang

y_radius = 5.5*units.Ang
x_radius = 5.0*units.Ang

#WIP DO NOT SET TRUE YET
mend_hexagons = False
# This is an offset for the center of the crack propagation,
# as the center is defined as the first broken bond, which is no quite what we want.

offset = [-1.0, 0.0, 0.0]
# ********** Setup calculator ************

renew_calculator = True 

# QMMM calculator from ASE, james kermode, and James's rescaled calculator
from atomistica import Tersoff, TersoffScr, Tersoff_PRB_39_5566_Si_C__Scr, Tersoff_PRB_39_5566_Si_C
mm_unscaled = Tersoff(**Tersoff_PRB_39_5566_Si_C)
#mm = mm_unscaled
mm = RescaledCalculator(TersoffScr(**Tersoff_PRB_39_5566_Si_C__Scr), 3.562, 528.624, 3.5656, 425.023)
#TersoffScr(**Tersoff_PRB_39_5566_Si_C__Scr)

##qm = SocketCalculator(castep_client)
qm = RescaledCalculator(TersoffScr(**Tersoff_PRB_39_5566_Si_C__Scr), 3.562, 528.624, 3.5656, 425.023)
qm_renew = RescaledCalculator(TersoffScr(**Tersoff_PRB_39_5566_Si_C__Scr), 3.562, 528.624, 3.5656, 425.023)
