#! /usr/bin/env python

# ======================================================================
# matscipy - Python materials science tools
# https://github.com/libAtoms/matscipy
#
# Copyright (2014) James Kermode, King's College London
#                  Lars Pastewka, Karlsruhe Institute of Technology
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ======================================================================

"""
Script to run classical molecular dynamics for a crack slab,
incrementing the load in small steps until fracture starts.
James Kermode <james.kermode@kcl.ac.uk>
August 2013

James Brixey <j.brixey@warwick.ac.uk>
December 2018
"""

from __future__ import print_function
import numpy as np

import ase.io
import ase.units as units
from ase import Atom
from ase.constraints import FixAtoms
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.netcdftrajectory import NetCDFTrajectory
#from ase.calculators.qmmm import ForceQMMM, RescaledCalculator
from qmmm import ForceQMMM, RescaledCalculator
from ase.neighborlist import neighbor_list
from matscipy.fracture_mechanics.crack import (get_strain,
                                               G_to_strain,
                                               strain_to_G,
                                               get_energy_release_rate,
                                               ConstantStrainRate,
                                               find_tip_broken_bonds,
                                               find_tip_stress_field)


from ase.optimize import FIRE, MDMin

from ase.geometry import get_distances
import ase.build

import sys
sys.path.insert(0, '.')
import diamond_thinstrip_crack_params as params

# ********** Read input file ************

print('Loading atoms from file "crack.xyz"')
atoms = ase.io.read('crack.xyz')

atoms *= (1, 1, 2)

orig_height = atoms.info['OrigHeight']
orig_crack_pos = atoms.info['CrackPos'].copy()

# ***** Setup constraints *******

top = atoms.positions[:, 1].max()
bottom = atoms.positions[:, 1].min()
left = atoms.positions[:, 0].min()
right = atoms.positions[:, 0].max()

# fix atoms in the top and bottom rows
fixed_mask = ((abs(atoms.positions[:, 1] - top) < 1.0) |
              (abs(atoms.positions[:, 1] - bottom) < 1.0))
fix_atoms = FixAtoms(mask=fixed_mask)
print('Fixed %d atoms\n' % fixed_mask.sum())

atoms.set_constraint(fix_atoms)
# Increase epsilon_yy applied to all atoms at constant strain rate

strain_atoms = ConstantStrainRate(orig_height,
                                  params.strain_rate*params.timestep)

atoms.set_constraint(fix_atoms)
atoms.set_calculator(params.mm_unscaled)

crack_origin = np.ndarray.tolist(find_tip_broken_bonds(atoms, 2.0))
#r = at0.get_distances(0, np.arange(1, len(at0)), mic=True)
atoms.set_calculator(params.mm)


# ********* Setup and run MD ***********

# Set the initial temperature to 2*simT: it will then equilibriate to
# simT, by the virial theorem
MaxwellBoltzmannDistribution(atoms, 2.0*params.sim_T)

# Initialise the dynamical system
dynamics = VelocityVerlet(atoms, params.timestep)

# Print some information every time step
def printstatus():
    if dynamics.nsteps == 1:
        print("""
State      step number    Temp/K     Strain      G/(J/m^2)  CrackPos/A D(CrackPos)/A
---------------------------------------------------------------------------------""")

    log_format = ('%(label)-4s%(time)12.1f%(temperature)12.6f'+
                  '%(strain)12.5f%(G)12.4f%(crack_pos_x)12.2f    (%(d_crack_pos_x)+5.2f)')

    atoms.info['label'] = 'D'                  # Label for the status line
    atoms.info['time'] = dynamics.get_time()/units.fs
    atoms.info['temperature'] = (atoms.get_kinetic_energy() /
                                 (1.5*units.kB*len(atoms)))
    atoms.info['strain'] = get_strain(atoms)
    atoms.info['G'] = get_energy_release_rate(atoms)/(units.J/units.m**2)

    atoms.set_calculator(param.mm_unscaled)
    crack_pos = find_tip_broken_bonds(atoms, 2.0)
    atoms.set_calculator(param.mm)
    atoms.info['crack_pos_x'] = crack_pos[0]
    atoms.info['d_crack_pos_x'] = crack_pos[0] - orig_crack_pos[0]

    print(log_format % atoms.info)


dynamics.attach(printstatus)

# Check if the crack has advanced enough and apply strain if it has not
def check_if_crack_advanced(atoms):
    atoms.set_calculator(params.mm_unscaled)
    crack_pos = find_tip_broken_bonds(atoms, 2.0)
    atoms.set_calculator(params.mm)

    # strain if crack has not advanced more than tip_move_tol
    if crack_pos[0] - orig_crack_pos[0] < params.tip_move_tol:
        strain_atoms.apply_strain(atoms)

dynamics.attach(check_if_crack_advanced, 1, atoms)

# Save frames to the trajectory every `traj_interval` time steps
#trajectory = NetCDFTrajectory(params.traj_file, mode='w')
#def write_frame(atoms):
#    trajectory.write(atoms)

# dynamics.attach(write_frame, diamond_thinstrip_crack_params.traj_interval, atoms)



E = 1281.2
nu = -0.001

# Start running!
#dynamics.run(diamond_thinstrip_crack_params.nsteps)

strain = G_to_strain(params.initial_G, E, nu, orig_height)

dump = open("crack_full_run.xyz", "w")
def write_frame(atoms=atoms):
    ase.io.write("crack_full_run.xyz", atoms, append=True)

atoms.new_array('qm_region', np.zeros(len(atoms)))
atoms.new_array('buffer_region', np.zeros(len(atoms)))

#traj = open("traj.xyz", "w")
#atoms_copy = atoms.copy()
#ase.io.write(traj, atoms_copy, format="xyz")

which_try = 1

for i in range(params.nsteps):

    #print(crack_origin)
    crack_structure = atoms.get_positions()
    cell = atoms.get_cell()
    pbc = atoms.get_pbc()
    print('pbc before', pbc)
    atoms.set_pbc([False,False,True])
    pbc = atoms.get_pbc()
    print('pbc after', pbc)

    if which_try==1:
        dummy_h = Atom(1, crack_origin)
        atoms.append(dummy_h)
        r = atoms.get_distances(len(atoms)-1, np.arange(0, len(atoms)-1), mic=True, vector=True)
	#print(r)
        
	del atoms[-1]
        r_prime = np.sqrt(r[:,0]**2+r[:,1]**2)
        mask = r_prime < params.quantum_region_size
	print('this is the mask')
	print('mask', mask)
#        atoms.set_mask() = mask
        print('QM atoms', mask.sum())
        atoms.arrays["qm_region"] = mask.astype(int)

        
    qmmm = ForceQMMM(atoms, mask, params.qm, params.mm, buffer_width=params.buffer_width, hydrogenate=False)
    atoms.set_calculator(qmmm)
    atoms.set_pbc([False, False, True])
    if qmmm.qm_buffer_mask is None:
        qmmm.initialize_qm_buffer_mask(atoms)
    mask1 = qmmm.qm_buffer_mask
    #np.set_printoptions(threshold=np.inf)
    print('this is mask 1')
    print(mask1)
    atoms.arrays["buffer_region"] = mask1.astype(int)
    #print(qmmm.qm_buffer_mask)
    #atoms.arrays["buffer_region"] = qmmm.qm_buffer_mask.astype(int)

    # ****** Apply initial strain ramp *****
    
    #if i == 1:
    #    atoms.positions[:, 1] *= ((1.0 + G_to_strain(9.0,  E, nu, orig_height)) / (1 + strain))

    new_strain = strain + params.strain_rate
    atoms.positions[:, 1] *= ((1.0 + new_strain)/(1 + strain))
    strain = new_strain

    print('Applied initial load: strain=%.4f, G=%.2f J/m^2' %
          (strain, strain_to_G(strain, E, nu,orig_height) / (units.J / units.m**2)))

    #f = atoms.get_forces()
    #atoms.set_array('qmmm_forces', f)
    #atoms_copy = atoms.copy()
    #ase.io.write(traj, atoms_copy, format="xyz")
    #traj.close()
    #1/0


    # ***** Relaxation of crack slab  *****

    # optionally, relax the slab, keeping top and bottom rows fixed
    if hasattr(params, 'relax_slab') and params.relax_slab:


        atoms.arrays["buffer_region"] = mask1.astype(int)
        #print(qmmm.qm_buffer_mask)
        #atoms.arrays["buffer_region"] = qmmm.qm_buffer_mask.astype(int)
        print('Relaxing slab...')
        #opt = LBFGSLineSearch(atoms)
        #opt.run(fmax=diamond_thinstrip_crack_params.relax_fmax_first)
        qmmm = ForceQMMM(atoms, mask, params.qm, params.mm, buffer_width=params.buffer_width, hydrogenate=True)
        atoms.set_calculator(qmmm)
        opt = FIRE(atoms)
        opt.attach(write_frame)
        opt.run(fmax=params.relax_fmax, steps = 10)
        #ase.io.write(all_opt_traj, atoms, format="xyz")
    #atoms_0 = atoms.copy()
    #atoms_0.set_calculator(params.mm)
    crack_pos = find_tip_broken_bonds(atoms, 2.0)
    if crack_pos[0] - orig_crack_pos[0] > 3:
        print('finished')
        break

    #mask2 = qmmm.qm_buffer_mask
    #atoms.arrays["buffer_region"] = mask2.astype(int)
    #atoms_copy = atoms.copy()
    #ase.io.write(traj, atoms_copy, format="xyz")

#traj.close()

