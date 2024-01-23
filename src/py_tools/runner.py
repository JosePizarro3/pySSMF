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

import logging
import os

from .input import ValidLatticeModels
from .schema import Model
from .parsing import MinimalWannier90Parser, ToyModels
from .hopping_pruning import Pruner
from .tb_hamiltonian import TBHamiltonian
from .visualization import plot_hopping_matrices, plot_band_structure


class Runner(ValidLatticeModels):
    def __init__(self, **kwargs):
        super().__init__()
        self.data = kwargs.get("data", {})
        # Initializing the model class
        self.model = Model()

    def parse_tb_model(self):
        lattice_model = self.data.get("tb_model", "")
        if "_hr.dat" in lattice_model:
            model_file = os.path.join(self.data.get("working_directory"), lattice_model)
            MinimalWannier90Parser().parse(model_file, self.model, self.logger)
        elif lattice_model in self._valid_lattice_models:
            onsite_energies = self.data.get("onsite_energies")
            hoppings = self.data.get("hoppings")
            ToyModels().parse(
                lattice_model, onsite_energies, hoppings, self.model, self.logger
            )
        else:
            self.logger.error(
                "Could not recognize the input tight-binding model. Please "
                "check your inputs."
            )

    def prune_hoppings(self):
        prune_threshold = self.data.get("prune_threshold")
        if prune_threshold:
            pruner = Pruner(self.model)
            pruner.prune_by_threshold(prune_threshold, self.logger)
            if self.data.get("plot_hoppings"):
                plot_hopping_matrices(pruner.hopping_matrix_norms / pruner.max_value)

    def calculate_band_structure(self):
        n_k_path = self.data.get("n_k_path", 90)
        tb_hamiltonian = TBHamiltonian(
            self.model, k_grid_type="bands", n_k_path=n_k_path
        )
        special_points = tb_hamiltonian.k_path.special_points
        kpoints = tb_hamiltonian.kpoints
        _, eigenvalues = tb_hamiltonian.diagonalize(kpoints)
        plot_band_structure(eigenvalues, tb_hamiltonian, special_points)

    def bz_diagonalization(self):
        k_grid = self.data.get("k_grid", [1, 1, 1])
        tb_hamiltonian = TBHamiltonian(self.model, k_grid_type="full_bz", k_grid=k_grid)
        kpoints = tb_hamiltonian.kpoints
        eigenvectors, eigenvalues = tb_hamiltonian.diagonalize(kpoints)

    def run(self):
        self.parse_tb_model()

        self.prune_hoppings()

        if self.data.get("plot_bands"):
            self.calculate_band_structure()

        self.bz_diagonalization()
