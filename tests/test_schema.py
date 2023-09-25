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


def test_system(example_system):
    '''Tests whether an instance of System is created and properly populated.'''
    assert example_system.n_atoms == 2
    assert example_system.labels == ['H', 'O']
    assert np.array_equal(
        example_system.positions.magnitude, np.array([[0.0, 0.0, 0.0], [0.1, 0.1, 0.1]]))
    assert np.array_equal(example_system.lattice_vectors.magnitude, np.identity(3))


def test_bravais_lattice(example_bravais_lattice):
    '''Tests whether an instance of BravaisLattice is created and properly populated.'''
    assert example_bravais_lattice.n_points == 3
    assert np.array_equal(
        example_bravais_lattice.points.magnitude,
        np.array([[0.0, 0.0, 0.0], [0.1, 0.1, 0.1], [0.2, 0.2, 0.2]]))
    assert example_bravais_lattice.system


def test_model(example_model):
    '''Tests whether an instance of Model is created and properly populated.'''
    assert example_model.n_orbitals == 4
    assert np.array_equal(example_model.degeneracy_factors, np.array([1, 2, 2]))
    assert np.array_equal(
        example_model.onsite_energies.magnitude, np.array([0.5, 0.6, 0.7, 0.8]))
    assert np.array_equal(example_model.hopping_matrix.magnitude, np.array([
        [[0.6, 0.5], [0.4, 0.3]],
        [[0.5, 0.4], [0.3, 0.2]],
        [[0.4, 0.3], [0.2, 0.1]]
    ]))
