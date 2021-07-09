#  Hybrid MD decision making package
#
#  Copyright (c) Tamas K. Stenczel 2021.
from ase.build import molecule
import ase.io
import numpy as np

cell_a = 15.0

# these should be centred at (0, 0, 0)
c60 = molecule("C60")
h2 = molecule("H2")

# move CoM to origin
c60.positions -= np.mean(c60.positions, axis=0)
h2.positions -= np.mean(h2.positions, axis=0)

# create a cell
atoms = h2 + c60
atoms.positions += cell_a / 2
atoms.cell = [cell_a, cell_a, cell_a]
atoms.set_pbc(True)  # dummy really

# write the cell file
ase.io.write("h2c60.cell", atoms)

# extra stuff in .cell
extra_in_cell = """
# calculation settings for gamma point only
kpoints_MP_grid 1 1 1
kpoints_MP_offset 0.0 0.0 0.0

# geometry optimisation specific
FIX_ALL_CELL: True
FIX_COM: True
"""

with open("h2c60.cell", "a") as file:
    file.write(extra_in_cell)
