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
import logging

from .schema import Model


class Pruner:
    def __init__(self, model: Model):
        """
        Initializes the `Pruner` object for the `Model` object.

        Args:
            model (Model): (Tight-binding) Model object to be pruned.
        """
        if not model:
            return
        self.model = model
        self.hopping_matrix_norms = None
        self.max_value = None
        self.update_norms_and_max_value()

    def update_norms_and_max_value(self):
        """
        Stores the hopping_matrix norms and their maximum value.
        """
        if self.model.hopping_matrix is not None:
            self.hopping_matrix_norms = np.abs(self.model.hopping_matrix.magnitude)
            self.max_value = np.max(self.hopping_matrix_norms)

    def prune_by_threshold(
        self, threshold_factor: float, logger: logging.Logger = None
    ):
        """
        Prune the model's hopping_matrix based on a threshold.

        Args:
            threshold_factor (float, optional): Percentage of the max_value to determine
                the pruning threshold. Defaults to 1% of the max_value (defined in input reading).
            logger (logging.Logger, optional): Logger object for debug messages. Defaults to None.
        """
        self.logger = logging.getLogger(__name__) if logger is None else logger
        # TODO improve this method. Right now it is assuming that hoppings are ordered from
        # larger to smaller values as far as the Bravais norm increases. This might give
        # problems for hoppings structures which are not so trivial.
        if not self.max_value and not self.hopping_matrix_norms:
            self.logger.warning(
                'Could not extract the hopping_matrix norms and their max_value.'
            )
            return
        threshold = threshold_factor * self.max_value

        matrix_sums = np.sum(self.hopping_matrix_norms, axis=(1, 2))
        n_orbitals = self.model.n_orbitals if self.model.n_orbitals else 1
        small_matrix_indices = np.where(
            matrix_sums < threshold * n_orbitals * n_orbitals
        )[0]
        if small_matrix_indices.size > 0:
            if small_matrix_indices[0] == 0 and small_matrix_indices[1] != 1:
                last_small_index = small_matrix_indices[1]
            else:
                last_small_index = small_matrix_indices[0]

            # Update model attributes
            try:
                self.model.hopping_matrix = self.model.hopping_matrix[:last_small_index]
                self.model.bravais_lattice.points = self.model.bravais_lattice.points[
                    :last_small_index
                ]
                self.model.bravais_lattice.n_points = last_small_index
                self.model.degeneracy_factors = self.model.degeneracy_factors[
                    :last_small_index
                ]
            except Exception:
                self.logger.warning(
                    'Could not update the model parameters after pruning.'
                )

        self.update_norms_and_max_value()
