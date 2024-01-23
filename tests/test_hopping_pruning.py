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
from src.hopping_pruning import Pruner


def test_pruner(example_model):
    """Tests whether the Pruner properly prunes a hopping_matrix quantity."""
    pruner = Pruner(example_model)

    # Initialization
    assert pruner.model == example_model
    assert pruner.max_value == 0.6

    # Pruning by input threshold
    pruner.prune_by_threshold(0.05)
    assert example_model.bravais_lattice.n_points == 3
    assert example_model.hopping_matrix.shape[0] == 3
    assert example_model.degeneracy_factors.shape[0] == 3
    assert np.array_equal(
        example_model.bravais_lattice.points[-1].magnitude, np.array([0.2, 0.2, 0.2])
    )

    pruner = Pruner(example_model)
    pruner.prune_by_threshold(0.15)
    assert example_model.bravais_lattice.n_points == 1
    assert example_model.hopping_matrix.shape[0] == 1
    assert example_model.degeneracy_factors.shape[0] == 1
    assert np.array_equal(
        example_model.bravais_lattice.points[-1].magnitude, np.array([0.0, 0.0, 0.0])
    )
