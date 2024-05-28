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

import pytest
import numpy as np

from src.schema import System, BravaisLattice, Model


def get_template_system():
    system = System()
    system.n_atoms = 2
    system.labels = ['H', 'O']
    system.positions = np.array([[0.0, 0.0, 0.0], [0.1, 0.1, 0.1]])
    system.lattice_vectors = np.identity(3)
    return system


def get_template_bravais_lattice():
    bravais_lattice = BravaisLattice()
    bravais_lattice.n_points = 3
    bravais_lattice.points = np.array(
        [[0.0, 0.0, 0.0], [0.1, 0.1, 0.1], [0.2, 0.2, 0.2]]
    )
    return bravais_lattice


def get_template_model():
    model = Model()
    model.n_orbitals = 4
    model.degeneracy_factors = np.array([1, 2, 2])
    model.onsite_energies = np.array([0.5, 0.6, 0.7, 0.8])
    model.hopping_matrix = np.array(
        [[[0.6, 0.5], [0.4, 0.3]], [[0.5, 0.4], [0.3, 0.2]], [[0.4, 0.3], [0.2, 0.1]]]
    )
    return model


@pytest.fixture
def example_system():
    return get_template_system()


@pytest.fixture
def example_bravais_lattice():
    system = get_template_system()
    bravais_lattice = get_template_bravais_lattice()
    bravais_lattice.system = system
    return bravais_lattice


@pytest.fixture
def example_model():
    system = get_template_system()
    bravais_lattice = get_template_bravais_lattice()
    bravais_lattice.system = system
    model = get_template_model()
    model.bravais_lattice = bravais_lattice
    return model
