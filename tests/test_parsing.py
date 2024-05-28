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
import os

from pyssmf.schema import Model
from pyssmf.parsing import MinimalWannier90Parser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


def test_wannier90_parser():
    """Tests the MinimalWannier90Parser with provided wannier90 files."""

    model = Model()
    filepath = os.path.join(
        os.path.dirname(__file__), 'data/wannier90/wannier90_hr.dat'
    )
    MinimalWannier90Parser().parse(filepath, model, None)

    # Model
    assert model.n_orbitals == 4
    assert model.degeneracy_factors.shape == (39,)
    assert np.array_equal(model.degeneracy_factors[:4], np.array([1, 1, 1, 1]))
    assert len(model.onsite_energies) == model.n_orbitals
    assert np.array_equal(
        model.onsite_energies.magnitude,
        np.array([-0.33423, -0.334384, -0.33423, -0.334384]),
    )
    assert model.hopping_matrix.shape == (39, model.n_orbitals, model.n_orbitals)
    assert model.hopping_matrix[2][0][0].magnitude == approx(-0.004462 + 0.000147j)

    # Bravais lattice
    bravais_lattice = model.bravais_lattice
    assert bravais_lattice.n_points == 39
    assert np.array_equal(
        bravais_lattice.points[2].magnitude, np.array([12.63, 0.0, 0.0])
    )

    # System
    system = bravais_lattice.system
    assert system.labels[:4] == ['Nb', 'Nb', 'Ta', 'Ta']
    assert np.array_equal(
        system.positions[0].magnitude, np.array([0.0, 7.29193, 22.06006])
    )
    lattice_vectors = np.array(
        [[12.63, 0.0, 0.0], [-6.315, 10.937901, 0.0], [0.0, 0.0, 25.26]]
    )
    assert np.array_equal(system.lattice_vectors.magnitude, lattice_vectors)
