from ase.io import read, write
from crack_matscipy_update import find_tip_broken_bonds

whole_set_of_atoms = read('tip_position_restricted_run_1.xyz', index = ':')


for i in range(0, len(whole_set_of_atoms)):
    new_atoms = read('tip_position_restricted_run_1.xyz', index = i)

    crack_origin = find_tip_broken_bonds(new_atoms, 2.25, 4)
    print(new_atoms)
    print(crack_origin)
    crack_origin_array = np.append(crack_origin_array, crack_origin)

    crack_origin = np.ndarray.tolist(find_tip_broken_bonds(atoms, 2.25, 4))
    #new_atoms.arrays["inner_section"] = fixed_mask.astype(int)
    write('tip_position_re_run.xyz', new_atoms, append=True)

print(crack_origin)

