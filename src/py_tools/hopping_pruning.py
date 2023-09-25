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


class Pruner:
    def __init__(self, model):
        self.model = model
        self.update_norms_and_max_value()

    def update_norms_and_max_value(self):
        """Stores the hopping_matrix norms and their maximum value.
        """
        self.hopping_matrix_norms = np.abs(self.model.hopping_matrix.magnitude)
        self.max_value = np.max(self.hopping_matrix_norms)

    def prune_by_threshold(self, threshold_factor: float = 0.05):
        """Prune the model's hopping_matrix based on a threshold.

        Args:
            threshold_factor (float, optional): Percentage of the max_value to determine
                the pruning threshold. Defaults to 5% of the max_value.
        """
        threshold = threshold_factor * self.max_value

        matrix_sums = np.sum(self.hopping_matrix_norms, axis=(1, 2))
        small_matrix_indices = np.where(matrix_sums < threshold * self.model.n_orbitals * self.model.n_orbitals)[0]
        if small_matrix_indices.size > 0:
            last_small_index = small_matrix_indices[0]

            # Update model attributes
            self.model.hopping_matrix = self.model.hopping_matrix[:last_small_index]
            self.model.bravais_lattice.points = self.model.bravais_lattice.points[:last_small_index]
            self.model.bravais_lattice.n_points = last_small_index
            self.model.degeneracy_factors = self.model.degeneracy_factors[:last_small_index]

        self.update_norms_and_max_value()
