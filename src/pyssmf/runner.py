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
import os
from scipy.stats import norm

from .input import ValidLatticeModels
from .schema import Model
from .parsing import MinimalWannier90Parser, ToyModels
from .hopping_pruning import Pruner
from .tb_hamiltonian import TBHamiltonian
from .visualization import plot_hopping_matrices, plot_band_structure, plot_dos


class Runner(ValidLatticeModels):
    """
    Class that runs the calculation. It reads the input data from the input file and
    runs the SSMF calculation.
    """

    def __init__(self, **kwargs):
        super().__init__()
        self.data = kwargs.get('data', {})
        # Initializing the model class
        self.model = Model()

    def parse_tb_model(self):
        """
        Parses the tight-binding model from the input file. It can be obtained from a Wannier90
        tight-binding calculation or a toy lattice model.
        """
        lattice_model = self.data.get('tb_model', '')
        if '_hr.dat' in lattice_model:
            model_file = os.path.join(self.data.get('working_directory'), lattice_model)
            MinimalWannier90Parser().parse(model_file, self.model, self.logger)
        elif lattice_model in self._valid_lattice_models:
            onsite_energies = self.data.get('onsite_energies')
            hoppings = self.data.get('hoppings')
            ToyModels().parse(
                lattice_model, onsite_energies, hoppings, self.model, self.logger
            )
        else:
            self.logger.error(
                'Could not recognize the input tight-binding model. Please '
                'check your inputs.'
            )
            return
        self.logger.info('Tight-binding model parsed successfully!')

    def prune_hoppings(self):
        """
        Prunes the hopping matrices by setting to zero all values below a certain `prune_threshold`.
        """
        prune_threshold = self.data.get('prune_threshold')
        if prune_threshold:
            pruner = Pruner(self.model)
            pruner.prune_by_threshold(prune_threshold, self.logger)
            if self.data.get('plot_hoppings'):
                plot_hopping_matrices(pruner.hopping_matrix_norms / pruner.max_value)
            self.logger.info('Hopping pruning finished!')

    def calculate_band_structure(self):
        """
        Calculates the band structure of the tight-binding model in a given `n_k_path`.
        """
        n_k_path = self.data.get('n_k_path', 90)
        tb_hamiltonian = TBHamiltonian(
            self.model, k_grid_type='bands', n_k_path=n_k_path
        )
        special_points = tb_hamiltonian.k_path.special_points
        kpoints = tb_hamiltonian.kpoints
        eigenvalues, _ = tb_hamiltonian.diagonalize(kpoints)
        plot_band_structure(eigenvalues, tb_hamiltonian, special_points)
        self.logger.info('Band structure calculation finished!')

    def gaussian_convolution(
        self,
        energies,
        orbital_dos_histogram,
        width: float = 0.1,
        delta_energy: float = 0.01,
    ):
        """
        Convolutes / smoothes the histogram data with a Gaussian distribution function as defined
        in scipy.stats.norm. The mesh of energies is also expanded (with delta_energy in eV)
        to resolve better the Gaussians.

        Args:
            energies (np.array): array containing the energies. Dimensions are (n_energies).
            orbital_dos_histogram (np.array): array containing the orbital DOS histogram. Dimensions
                are (n_energies, n_orbitals).
            width (float): standard deviation or width of the Gaussian distribution. Defaults to
                0.1 eV.
            delta_energy (float): the spacing of the new energies mesh. Defaults to 0.01 eV.

        Returns:
            new_energies, convoluted_data: returns the new X and Y data for the convoluted data.
        """
        energy_min = np.min(energies) - 2 * width
        energy_max = np.max(energies) + 2 * width
        new_energies = np.arange(energy_min, energy_max, delta_energy)

        n_orbitals = orbital_dos_histogram.shape[1]

        gaussian = norm(loc=0, scale=width)

        convoluted_data = np.zeros((len(new_energies), n_orbitals))
        for i, new_energy in enumerate(new_energies):
            for j in range(n_orbitals):
                convoluted_data[i, j] = np.sum(
                    orbital_dos_histogram[:, j] * gaussian.pdf(new_energy - energies)
                )
        return new_energies, np.transpose(convoluted_data)

    def calculate_dos(
        self,
        eigenvalues,
        eigenvectors,
        bins: int = 100,
        width: float = 0.1,
        delta_energy: float = 0.01,
    ):
        """
        Calculates the orbital and total DOS from the eigenvalues and eigenvectors of the
        tight-binding model.

        Args:
            eigenvalues (np.ndarray): the eigenvalues of the tight-binding model. Dimensions are
                (n_kpoints, n_orbitals).
            eigenvectors (np.ndarray): the eigenvectors of the tight-binding model. Dimensions are
                (n_kpoints, n_orbitals, n_orbitals).
            bins (int, optional): the number of bins for the histogram. Defaults to 100.
            width (float, optional): the width of the Gaussian distribution function. Defaults to
                0.1 eV.
            delta_energy (float, optional): the spacing of the new energies mesh. Defaults to
                0.01 eV.

        Returns:
            energies, orbital_dos, total_dos: the energies mesh, the orbital DOS, and the total DOS.
        """
        # We create the orbital DOS histogram and append all orbital contributions together
        orbital_dos_histogram = []
        for orbital in range(self.model.n_orbitals):
            orbital_dos_contribution_histogram, bin_edges = np.histogram(
                eigenvalues,
                bins=bins,
                density=True,
                weights=np.abs(eigenvectors[:, orbital]) ** 2,
            )
            energies = (bin_edges[:-1] + bin_edges[1:]) / 2
            orbital_dos_histogram.append(orbital_dos_contribution_histogram)
        if not orbital_dos_histogram:
            self.logger.warning(
                'Problem obtaining the orbital DOS histogram. Cannot resolve DOS.'
            )
        # We convolute the histogram to obtain a smoother orbital DOS
        energies, orbital_dos = self.gaussian_convolution(
            energies, np.array(orbital_dos_histogram).T, width, delta_energy
        )
        # We sum all orbital contributions to obtain the total DOS
        total_dos = np.sum(orbital_dos, axis=0)
        return energies, orbital_dos, total_dos

    def bz_diagonalization(self):
        """
        Diagonalizes the tight-binding model in the full Brillouin zone and returns its
        eigenvalues and eigenvectors.
        """
        k_grid = self.data.get('k_grid', [1, 1, 1])
        tb_hamiltonian = TBHamiltonian(self.model, k_grid_type='full_bz', k_grid=k_grid)
        kpoints = tb_hamiltonian.kpoints
        eigenvalues, eigenvectors = tb_hamiltonian.diagonalize(kpoints)

        # Calculating and plotting DOS
        if self.data.get('dos'):
            bins = int(np.linalg.norm(k_grid))
            width = self.data.get('dos_gaussian_width')
            delta_energy = self.data.get('dos_delta_energy')
            energies, orbital_dos, total_dos = self.calculate_dos(
                eigenvalues, eigenvectors, bins, width, delta_energy
            )
            plot_dos(energies, orbital_dos, total_dos)
            self.logger.info('DOS calculation finished!')

        self.logger.info('BZ diagonalization calculation finished!')
        return eigenvalues, eigenvectors

    def run(self):
        self.parse_tb_model()

        self.prune_hoppings()

        if self.data.get('plot_bands'):
            self.calculate_band_structure()

        self.bz_diagonalization()
