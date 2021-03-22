import uuid
import nglview

import numpy as np
import matplotlib.pyplot as plt

from ipywidgets import Dropdown, FloatSlider, VBox, Output


@nglview.register_backend('ase')
class MyASEStructure(nglview.Structure):
    def __init__(self, atoms, bfactor=[], occupancy=[]):
        # super(MyASEStructure, self).__init__()
        self.ext = 'pdb'
        self.params = {}
        self._atoms = atoms

        self.bfactor = bfactor  # [min, max]
        self.occupancy = occupancy  # [0, 1]

        self.id = str(uuid.uuid4())

    def get_structure_string(self):
        """Example PDB file format:
        CRYST1   16.980   62.517  124.864  90.00  90.00  90.00 P 1
        MODEL     1
        ATOM      0   Fe MOL     1      15.431  60.277   6.801  1.00  0.00          FE
        ATOM      1   Fe MOL     1       1.273   3.392  93.940  1.00  0.00          FE
        """

        data = ""

        if self._atoms.get_pbc().any():
            cellpar = self._atoms.get_cell_lengths_and_angles()

            str_format = 'CRYST1' + '{:9.3f}' * 3 + '{:7.2f}' * 3 + ' P 1\n'
            data += str_format.format(*cellpar.tolist())

        data += 'MODEL     1\n'

        str_format = 'ATOM  {:5d} {:>4s} MOL     1    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:2s}\n'
        for index, atom in enumerate(self._atoms):
            data += str_format.format(
                index,
                atom.symbol,
                atom.position[0].tolist(),
                atom.position[1].tolist(),
                atom.position[2].tolist(),
                self.occupancy[index] if index <= len(self.occupancy) - 1 else 1.0,
                self.bfactor[index] if index <= len(self.bfactor) - 1 else 1.0,
                atom.symbol.upper()
            )

        data += 'ENDMDL\n'

        return data


def ViewStructure(atoms, repetition=(1, 1, 1)):
    # visualisation
    view = nglview.NGLWidget(gui=False)

    view.stage.set_parameters(clip_dist=0)
    view.add_structure(MyASEStructure(atoms))
    view.add_unitcell()

    view.add_structure(MyASEStructure(atoms.repeat(repetition)))

    return view


class AtomViewer(object):
    def __init__(self, atoms, data=[], xsize=1000, ysize=500):
        self.view = self._init_nglview(atoms, data, xsize, ysize)

        self.widgets = {
            'radius': FloatSlider(
                value=0.8, min=0.0, max=1.5, step=0.01,
                description='Ball size'
            ),
            'color_scheme': Dropdown(description='Solor scheme:'),
            'colorbar': Output()
        }
        self.show_colorbar(data)

        self.widgets['radius'].observe(self._update_repr)

        self.gui = VBox([
            self.view,
            self.widgets['colorbar'],
            self.widgets['radius']])

    def _update_repr(self, chg=None):
        self.view.update_spacefill(
            radiusType='radius',
            radius=self.widgets['radius'].value
        )

    def show_colorbar(self, data):
        with self.widgets['colorbar']:
            # Have colormaps separated into categories:
            # http://matplotlib.org/examples/color/colormaps_reference.html
            cmap = 'rainbow'

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 2))
            img = ax1.imshow([[min(data), max(data)]], aspect='auto', cmap=plt.get_cmap(cmap))
            ax1.remove()
            cbar = fig.colorbar(img, cax=ax2, orientation='horizontal')

            plt.show()

    @staticmethod
    def _init_nglview(atoms, data, xsize, ysize):
        view = nglview.NGLWidget(gui=False)
        view._remote_call(
            'setSize',
            target='Widget',
            args=[
                '{:d}px'.format(xsize),
                '{:d}px'.format(ysize)
            ]
        )

        data = np.max(data) - data

        structure = MyASEStructure(atoms, bfactor=data)
        view.add_structure(structure)

        view.clear_representations()
        view.add_unitcell()

        view.add_spacefill(
            # radiusType='radius',
            # radius=1.0,
            color_scheme='bfactor',
            color_scale='rainbow'
        )
        view.update_spacefill(
            radiusType='radius',
            radius=1.0
        )

        # update camera type
        view.control.spin([1, 0, 0], np.pi / 2)
        view.control.spin([0, 0, 1], np.pi / 2)
        view.camera = 'orthographic'
        view.center()

        return view
