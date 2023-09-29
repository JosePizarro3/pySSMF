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

from .schema import Model
from .parsing import MinimalWannier90Parser, ToyModels
from .hopping_pruning import Pruner
from .tb_hamiltonian import TBHamiltonian
from .visualization import plot_hopping_matrices, plot_band_structure


class Runner:
    _valid_lattice_models = ['linear', 'square', 'honeycomb', 'triangular']

    def __init__(self, data):
        self.model = Model()
        self.data = data

    def parse_tb_model(self):
        lattice_model = self.data.get('model', '')
        if '_hr.dat' in lattice_model:
            model_file = os.path.join(self.data.get('working_directory'), lattice_model)
            MinimalWannier90Parser().parse(model_file, self.model, self.logger)
        elif lattice_model in self._valid_lattice_models:
            onsite_energies = self.data.get('onsite_energies')
            hoppings = self.data.get('hoppings')
            ToyModels().parse(lattice_model, onsite_energies, hoppings, self.model, self.logger)

    def prune_hoppings(self, plot_hoppings):
        prune_threshold = self.data.get('prune_threshold')
        if prune_threshold:
            pruner = Pruner(self.model)
            pruner.prune_by_threshold(prune_threshold, self.logger)
            if plot_hoppings:
                plot_hopping_matrices(pruner.hopping_matrix_norms / pruner.max_value)

    def calculate_band_structure(self):
        tb_hamiltonian = TBHamiltonian(self.model, k_grid_type='bands')
        special_points = tb_hamiltonian.k_path.special_points
        _, eigenvalues = tb_hamiltonian.diagonalize(tb_hamiltonian.kpoints)
        plot_band_structure(eigenvalues, tb_hamiltonian, special_points)

    def run(self, plot_hoppings: bool = False, plot_bands: bool = False, logger: logging.Logger = None):
        self.logger = logging.getLogger(__name__) if logger is None else logger

        self.parse_tb_model()

        self.prune_hoppings(plot_hoppings)

        if plot_bands:
            self.calculate_band_structure()