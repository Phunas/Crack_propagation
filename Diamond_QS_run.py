#y - Python materials science tools
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
import sys
sys.path.insert(0, '.')

import numpy as np

import ase.io
import ase.units as units
from ase import Atom
from ase.constraints import FixAtoms
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.netcdftrajectory import NetCDFTrajectory
from ase.neighborlist import neighbor_list
from ase.optimize import FIRE, MDMin
from ase.geometry import get_distances
import ase.build


from qmmm import ForceQMMM, RescaledCalculator
from crack_matscipy_update import (get_strain,
                                   G_to_strain,
                                   strain_to_G,
                                   get_energy_release_rate,
                                   ConstantStrainRate,
                                   find_tip_broken_bonds,
                                   find_tip_stress_field)

import params

# ********** Force Printing *************

def flush():
    sys.stdout.flush()


# ********** Read input file ************

print('Loading atoms from file "thermalisation_faster.xyz"')
atoms = ase.io.read('crack.xyz')

if params.defect:
    atoms[params.defect_location].number = params.defect_number

if params.vacancy:
    del(atoms[params.vacancy_location].number)
    


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

strain = get_strain(atoms)
print(strain)
atoms.set_constraint(fix_atoms)

# Finding the Crack Tip

crack_origin = np.ndarray.tolist(find_tip_broken_bonds(atoms, 2.25, 4))
print('crack tip: ', crack_origin)


#This is the initial strain set
#strain = G_to_strain(params.initial_G, params.E, params.nu, orig_height)

#These save the frame
dump = open("crack_run.xyz", "w")
def write_frame_2(atoms=atoms):
    ase.io.write("crack_run.xyz", atoms, append=True)



crack_structure = atoms.get_positions()
cell = atoms.get_cell()
atoms.set_pbc([False,False,True])

new_strain = G_to_strain(strain_to_G(strain, params.E, params.nu, orig_height) + 8.0*units.J/units.m**2,
                         params.E, params.nu, orig_height)
atoms.positions[:, 1] *= ((1.0 + new_strain)/(1 + strain))
strain = new_strain
new_strain = strain



which_try=1
if which_try==1:
    dummy_h = Atom(1, np.subtract(crack_origin,params.offset))
    atoms.append(dummy_h)
    r = atoms.get_distances(len(atoms)-1, np.arange(0, len(atoms)-1), mic=True, vector=True)
    del atoms[-1]
if which_try==2:
    r = ase.geometry.get_distances([crack_origin], [atoms.get_positiions()] ,cell = atoms.get_cell() ,pbc = atoms.get_pbc())
r_prime_elipse = r[:,0]**2/params.x_radius**2 + r[:,1]**2/params.y_radius**2
mask = r_prime_elipse < 1
print('this is the mask')
print('mask', mask)
print('QM atoms', mask.sum())
atoms.arrays["qm_region"] = mask.astype(int)


qmmm = ForceQMMM(atoms, mask, params.qm_renew, params.mm, buffer_width=params.buffer_width,
                     hydrogenate=True, save_clusters=True, fixed_buffer=False, atom_of_interest=crack_origin)

atoms.set_calculator(qmmm)
if qmmm.qm_buffer_mask is None:
    qmmm.initialize_qm_buffer_mask(atoms)

mask1 = qmmm.qm_buffer_mask
atoms.arrays["buffer_region"] = mask1.astype(int)


# Increase epsilon_yy applied to all atoms at constant strain rate


atoms.set_constraint(fix_atoms)


crack_origin = np.ndarray.tolist(find_tip_broken_bonds(atoms, 2.25, 4))
#r = at0.get_distances(0, np.arange(1, len(at0)), mic=True)

atoms.set_calculator(qmmm)
crack_origin_vec = []
def update_qm_region():

    which_try=1
    if which_try==1:
        crack_origin = np.ndarray.tolist(find_tip_broken_bonds(atoms, 2.25, 4))
        crack_origin_vec.append(crack_origin)
        if len(crack_origin) >=5:
            c = len(crack_origin_vec)
            crack_origin = (crack_origin_vec[c] + crack_origin_vec[c-1] + crack_origin_vec[c-2] + crack_origin_vec[c-3])/4
        dummy_h = Atom(1, np.subtract(crack_origin,params.offset))
        atoms.append(dummy_h)
        r = atoms.get_distances(len(atoms)-1, np.arange(0, len(atoms)-1), mic=True, vector=True)
        del atoms[-1]
    if which_try==2:
        crack_origin = np.ndarray.tolist(find_tip_broken_bonds(atoms, 2.25, 4))
        r = ase.geometry.get_distances([crack_origin], [atoms.get_positiions()] ,cell = atoms.get_cell() ,pbc = atoms.get_pbc())
    r_prime_elipse_in = r[:,0]**2/params.x_radius**2 + r[:,1]**2/params.y_radius**2
    r_prime_elipse_out = r[:,0]**2/(params.x_radius+1)**2 + r[:,1]**2/(params.y_radius+1)**2
 
    if dynamics.nsteps == 0:
        mask_p = (r_prime_elipse_in < 1)
    else:
        mask_p = dynamics.mask.copy()
    
    mask = (r_prime_elipse_in < 1) | (r_prime_elipse_out < 1 & mask_p)
    #
    #atoms.calc.qm = params.qm_reknew or qm
    # even better atoms.calc.qm.params.reuse = True/False



    qmmm = ForceQMMM(atoms, mask, params.qm, params.mm, buffer_width=params.buffer_width,
               hydrogenate=True, save_clusters=True, fixed_buffer=False, atom_of_interest=crack_origin)


    if dynamics.nsteps > 0:
        if params.renew_calculator == True:
            if not qmmm.is_compatible(atoms):
                print('default')
                qmmm.qm_calc.param.reuse = 'default'

    atoms_prev = atoms.copy()
    dynamics.mask = mask.copy()

    print('this is the mask')
    print('mask', mask)
    print('QM atoms', mask.sum())
    qmmm.qm_selection_mask = mask
    qmmm.initialize_qm_buffer_mask(atoms) 
    dynamics.cell = qmmm.cell
    print('this is the cell for the qm region maybe', qmmm.cell)

    mask1 = qmmm.qm_buffer_mask

    atoms.arrays["qm_region"] = mask.astype(int)
    atoms.arrays["buffer_region"] = mask1.astype(int)
    
    atoms.set_calculator(qmmm)

# ********* Setup and run MD ***********

# Set the initial temperature to 2*simT: it will then equilibriate to
# simT, by the virial theorem
#MaxwellBoltzmannDistribution(atoms, 600 * units.kB)


# Initialise the dynamical system
dynamics = VelocityVerlet(atoms, params.timestep)

# Print some information every time step
def printstatus():
    if dynamics.nsteps == 0:
        print("""
State      Time/fs    Temp/K     Strain      G/(J/m^2)  CrackPos/A D(CrackPos)/A
---------------------------------------------------------------------------------""")

    log_format = ('%(label)-4s%(time)12.1f%(temperature)12.6f'+
                  '%(strain)12.5f%(G)12.4f%(crack_pos_x)12.2f    (%(d_crack_pos_x)+5.2f)')

    atoms.info['label'] = 'D'                  # Label for the status line
    atoms.info['time'] = dynamics.get_time()/units.fs
    atoms.info['temperature'] = (atoms.get_kinetic_energy() /
                                 (1.5*units.kB*len(atoms)))
    atoms.info['strain'] = get_strain(atoms)
    atoms.info['G'] = get_energy_release_rate(atoms)/(units.J/units.m**2)

    crack_pos = find_tip_broken_bonds(atoms, 2.25, 4) 
    atoms.info['crack_pos_x'] = crack_pos[0]
    atoms.info['d_crack_pos_x'] = crack_pos[0] - orig_crack_pos[0]

    print(log_format % atoms.info)

dynamics.attach(printstatus)

# Check if the crack has advanced enough and apply strain if it has not
def check_if_crack_advanced(atoms):
    crack_pos = find_tip_broken_bonds(atoms, 2.35, 4)

    # strain if crack has not advanced more than tip_move_tol
    if crack_pos[0] - orig_crack_pos[0] < params.tip_move_tol:
        print('crack has advanced', crack_pos[0] - orig_crack_pos[0])
#        strain_atoms.apply_strain(atoms)

dynamics.attach(check_if_crack_advanced, 1, atoms)

# Save frames to the trajectory every `traj_interval` time steps
dynamics.attach(flush)
dynamics.attach(update_qm_region)
dynamics.attach(write_frame_2)

# Start running!
dynamics.run(params.nsteps)
