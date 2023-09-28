#
# Copyright: Dr. José M. Pizarro.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import numpy as np
import enum
from typing import Tuple

import ase
from ase.spacegroup import get_spacegroup, spacegroup
from ase.dft.kpoints import monkhorst_pack, BandPath

from nomad.atomutils import Formula
from nomad.units import ureg
from .schema import Model


class KSampling:
    def __init__(self, model: Model, k_grid: list):
        """Initializes the `KSampling` object for the `Model` object and the k_grid list to
        generate the `ase.Atoms` cell object and thus the: spacegroup, k_path, and k_mesh
        properties.

        Args:
            model (Model): Input tight-binding model class.
            k_grid (list): list of k_grid for generating the Monkhorst-Pack mesh, [x, y, z].
                Default (defined in `TBHamiltonian` is [1, 1, 1])
        """
        if not model:
            return
        self.model = model
        self.system = model.bravais_lattice.system
        self.atoms = self.set_ase_atoms()
        self.k_grid = k_grid

    def set_ase_atoms(self) -> ase.Atoms:
        """Sets the `ase.Atoms` cell object. It also calculates the reciprocal_lattice_vectors
        and the Hill formula of the system.

        Returns:
            ase.Atoms: A class of ase.Atoms filled up with the system information.
        """
        atom_labels = self.system.labels
        pbc = self.system.periodic
        lattice_vectors = self.system.lattice_vectors.magnitude
        positions = self.system.positions.magnitude
        atoms = ase.Atoms(
            symbols=atom_labels,
            pbc=pbc,
            cell=lattice_vectors,
            positions=positions
        )
        reciprocal_lattice_vectors = 2 * np.pi * atoms.get_reciprocal_cell()
        self.system.reciprocal_lattice_vectors = reciprocal_lattice_vectors / ureg.angstrom

        try:
            formula = Formula(atoms.get_chemical_formula())
            self.model.bravais_lattice.formula_hill = formula.format('hill')
        except Exception:
            pass

        return atoms

    @property
    def spacegroup(self) -> spacegroup.Spacegroup:
        """Returns the spacegroup of the system based on the `ase.Atoms` cell object.

        Returns:
            spacegroup.Spacegroup: A spacegroup object representing the crystallographic spacegroup.
        """
        return get_spacegroup(self.atoms)

    @property
    def k_path(self) -> BandPath:
        """Returns the k-path for the system based on the `ase.Atoms` cell object.

        Returns:
            BandPath: A BandPath object representing the k-path.
        """
        lattice = self.atoms.cell.get_bravais_lattice()
        special_points = lattice.get_special_points()
        k_points = [list(value) for value in special_points.values()]
        return self.atoms.cell.bandpath(k_points, npoints=90)

    @property
    def k_mesh(self) -> np.ndarray:
        """Returns the Monkhorst-Pack k-mesh for the system.

        Returns:
            np.ndarray: A np.array representing the Monkhorst-Pack k-mesh.
        """
        return monkhorst_pack(self.k_grid) @ self.system.reciprocal_lattice_vectors.magnitude


class TBHamiltonian(KSampling):
    _valid_k_grid_types = ['bands', 'full_bz']

    def __init__(self, model: Model, k_grid_type: str, k_grid: list = [1, 1, 1]):
        """Initializes the `TBHamiltonian` object for the `Model` object, k_grid_type, and optional k_grid
        for generating kpoints. It also calculates the Hamiltonian matrix dimensions
        (n_k_points, n_orbitals, n_orbitals), and the number of R Bravais lattice vectors,
        n_r_points.

        Args:
            model (Model): Input tight-binding model class.
            k_grid_type (str): Enum with the type of k-point grid ('bands' or 'full_bz').
            k_grid (list, optional): List of k_grid for generating the Monkhorst-Pack mesh
                for the 'full_bz' calculation. Default is [1, 1, 1].
        """
        if k_grid_type not in self._valid_k_grid_types:
            raise ValueError("Invalid k_grid_type. Please, use 'bands' or 'full_bz'.")
        super().__init__(model, k_grid)
        self.k_grid_type = k_grid_type
        self.kpoints = np.empty((0, 3))  # initializing for mypy
        if k_grid_type == 'bands':
            self.kpoints = self.k_path.cartesian_kpts()
        elif k_grid_type == 'full_bz':
            self.kpoints = self.k_mesh
        self.n_orbitals = self.model.n_orbitals
        self.n_k_points = len(self.kpoints)
        self.n_r_points = self.model.bravais_lattice.n_points

    def __repr__(self) -> str:
        cls_name = self.__class__.__name__
        args = [
            f'n_orbitals={self.n_orbitals}',
            f'n_k_points={self.n_k_points}',
            f'n_r_points={self.n_r_points}',
            f'k_grid_type={self.k_grid_type}'
        ]
        if self.spacegroup:
            args.append(f'spacegroup.no={self.spacegroup.no}')
        return f"{cls_name}({', '.join(filter(None, args))})"

    def hamiltonian(self, kpoints: np.ndarray) -> np.ndarray:
        """Returns the Hamiltonian matrix for given k-points.

        Args:
            kpoints (np.ndarray): Array of k-points at which to calculate the Hamiltonian.

        Returns:
            np.ndarray: The Hamiltonian matrix for the specified k-points.
        """
        n_orbitals = self.model.n_orbitals
        n_rpoints = self.model.bravais_lattice.n_points
        onsite_energies = np.array([self.model.onsite_energies.magnitude])  # use to define same shape as hopping_matrix
        hamiltonian = np.zeros((len(kpoints), n_orbitals, n_orbitals), dtype=complex)
        for nr in range(n_rpoints):
            r_vector = self.model.bravais_lattice.points.magnitude[nr]
            exp_factor = np.exp(- 1j * np.pi * np.dot(kpoints, r_vector))
            hop = self.model.hopping_matrix.magnitude[nr] if nr != 0 else onsite_energies
            deg_factor = self.model.degeneracy_factors[nr]
            hamiltonian += hop[np.newaxis, :, :] * exp_factor[:, np.newaxis, np.newaxis] / deg_factor
        return hamiltonian

    def diagonalize(self, kpoints: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Diagonalizes the Hamiltonian matrix for given k-points and returns its eigenvectors
        and eigenvalues.

        Args:
            kpoints (np.ndarray): Array of k-points at which to diagonalize the Hamiltonian.

        Returns:
            Tuple[np.ndarray, np.ndarray]: A tuple containing eigenvectors and eigenvalues.
        """
        eigenvalues, eigenvectors = np.linalg.eigh(self.hamiltonian(kpoints))
        return eigenvectors, eigenvalues
