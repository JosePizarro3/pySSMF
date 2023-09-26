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
import subprocess

import ase
from ase.spacegroup import get_spacegroup

from nomad.atomutils import Formula
from nomad.units import ureg


class KSampling:
    def __init__(self, model):
        if not model:
            return
        self.model = model
        self.system = model.bravais_lattice.system
        self.atoms = self.set_ase_atoms()

    def set_ase_atoms(self):
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
        reciprocal_lattice_vectors = atoms.get_reciprocal_cell()
        self.system.reciprocal_lattice_vectors = reciprocal_lattice_vectors / ureg.angstrom

        try:
            formula = Formula(atoms.get_chemical_formula())
            self.model.bravais_lattice.formula_hill = formula.format('hill')
        except Exception:
            pass

        return atoms

    @property
    def spacegroup(self):
        return get_spacegroup(self.atoms)

    @property
    def k_path(self):
        lattice = self.atoms.cell.get_bravais_lattice()
        special_points = lattice.get_special_points()
        k_points = [list(value) for value in special_points.values()]
        return self.atoms.cell.bandpath(k_points, npoints=90)


class TBHamiltonian(KSampling):
    def __init__(self, model):
        super().__init__(model)
        self.kpoints = self.k_path.cartesian_kpts()
        self.n_orbitals = self.model.n_orbitals
        self.n_k_points = len(self.kpoints)
        self.n_r_points = self.model.bravais_lattice.n_points
    
    def __repr__(self):
        cls_name = self.__class__.__name__
        args = [
            f'n_orbitals={self.n_orbitals}',
            f'n_k_points={self.n_k_points}',
            f'n_r_points={self.n_r_points}',
            f'k_path=<{self.k_path.path}>',
            f'spacegroup=<{self.spacegroup.no}>'
        ]
        return f"{cls_name}({', '.join(filter(None, args))})"

    def hamiltonian(self, kpoints):
        n_orbitals = self.model.n_orbitals
        n_rpoints = self.model.bravais_lattice.n_points
        onsite_energies = np.array([self.model.onsite_energies.magnitude])  # use to define same shape as hopping_matrix
        hamiltonian = np.zeros((len(kpoints), n_orbitals, n_orbitals), dtype=complex)
        for nr in range(n_rpoints):
            r_vector = self.model.bravais_lattice.points.magnitude[nr]
            exp_factor = np.exp(- 1j * np.pi * np.dot(kpoints, r_vector))
            hop = self.model.hopping_matrix.magnitude[nr] if nr != 0 else onsite_energies
            hamiltonian += hop[np.newaxis, :, :] * exp_factor[:, np.newaxis, np.newaxis]
        return hamiltonian
    
    def diagonalize(self, kpoints):
        eigenvalues, eigenvectors = np.linalg.eigh(self.hamiltonian(kpoints))
        return eigenvectors, eigenvalues
