from atomistica import TightBinding
from ase.lattice import bulk
import ase.units as units
import ase.build
#from hotbit import Hotbit
from qmmm import RescaledCalculator

# ********** Bulk unit cell ************

# 8-atom diamond cubic unit cell for diamond
a0_cast = 3.5694 # guess at lattice constant for carbon - we will minimize
a0_tersoff = 3.5
#cryst = bulk('C', 'diamond', a=a0, cubic=True)
cryst = ase.build.bulk('C', 'diamond', a=a0_cast, cubic=True)
surf_ny = 4 # number of cells needed to get accurate surface energy

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
depth = 1			     # The number of layers in the slab
vacuum = 20.0*units.Ang             # Amount of vacuum around slab
crack_seed_length = 50.0*units.Ang   # Length of seed crack
strain_ramp_length = 10.0*units.Ang  # Distance over which strain is ramped up
initial_G = 10.0*(units.J/units.m**2) # Initial energy flow to crack tip

relax_bulk = True                     # If True, relax initial bulk cell
bulk_fmax  = 1e-6*units.eV/units.Ang  # Max force for bulk, C_ij and surface energy

relax_slab = True                     # If True, relax notched slab with calculator
relax_fmax = 0.01*units.eV/units.Ang # Maximum force criteria for relaxation

# ******* Molecular dynamics parameters ***********

sim_T = 300.0*units.kB           # Simulation temperature
nsteps = 1                       # Total number of timesteps to run for
timestep = 1.0*units.fs          # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang      # Amount by which potential cutoff is increased
                                 # for neighbour calculations
tip_move_tol = 10.0              # Distance tip has to move before crack
                                 # is taken to be running
strain_rate = 1e-3#*(1/units.fs)  # Strain rate
traj_file = 'traj.nc'            # Trajectory output file (NetCDF format)
traj_interval = 10               # Number of time steps between
                                 # writing output frames

# ******** Multiscale parameters *********

buffer_width = 3.0*units.Ang

quantum_region_size = 3.0*units.Ang

# ********** Setup calculator ************

# Stillinger-Weber (SW) classical interatomic potential, from QUIP
#from quippy import Potential
#calc = Potential('IP SW', 'params.xml')

# Screened Kumagai potential, from Atomistica
#import atomistica
#calc = atomistica.KumagaiScr()

# QMMM calculator from ASE, james kermode, and James's rescaled calculator
from atomistica import Tersoff, TersoffScr, Tersoff_PRB_39_5566_Si_C__Scr, Tersoff_PRB_39_5566_Si_C
mm_unscaled = Tersoff(**Tersoff_PRB_39_5566_Si_C)
#mm = mm_unscaled
mm = RescaledCalculator(Tersoff(**Tersoff_PRB_39_5566_Si_C), 3.562, 528.624, 3.5656, 425.023) 
from ase.calculators.dftb import Dftb

#qm = Hotbit(SCC=False, gamma_cut=5.0, verbose_SCC=True, width=0.05, txt='Vacancy.cal', kpts=(1,1,6), elements={'C':'pbc-0-3/C-C.skf', 'H' : 'pbc-0-3/H-H.skf'},
#               tables = {'CC':'pbc-0-3/C-C.skf', 'CH':'pbc-0-3/C-H.skf', 'HH' : 'pbc-0-3/H-H.skf'})

#from atomistica import TightBinding

qm = TightBinding(width=0.05, database_folder='./pbc-0-3')#, elements={'C':'pbc-0-3/C-C.skf', 'H' : 'pbc-0-3/H-H.skf'})#,

from ase.calculators.lj import LennardJones
sigma = 2.9035*(2.**(-1./6))

#qm = LennardJones(sigma=sigma, epsilon=0.05)
#from ase.build import bulk
#c = bulk("C")
#c.set_calculator(qm)
#print c.get_forces()
