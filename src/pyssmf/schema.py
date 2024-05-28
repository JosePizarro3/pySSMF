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

from nomad.metainfo import MSection, Quantity, Section, SubSection


class System(MSection):
    """
    Section containing the atomic system metadata.
    """

    m_def = Section(validate=False)

    n_atoms = Quantity(
        type=np.int32,
        description="""
        The total number of atoms in the system.
        """,
    )

    labels = Quantity(
        type=str,
        shape=["n_atoms"],
        description="""
        List containing the labels of the atoms. In the usual case, these correspond to
        the chemical symbols of the atoms. One can also append an index if there is a
        need to distinguish between species with the same symbol, e.g., atoms of the
        same species assigned to different atom-centered basis sets or pseudo-potentials,
        or simply atoms in different locations in the structure such as those in the bulk
        and on the surface. In the case where a species is not an atom, and therefore
        cannot be representated by a chemical symbol, the label can simply be the name of
        the particles.
        """,
    )

    positions = Quantity(
        type=np.float64,
        shape=["n_atoms", 3],
        unit="angstrom",
        description="""
        Positions of all the species, in cartesian coordinates. This metadata defines a
        configuration and is therefore required. For alloys where concentrations of
        species are given for each site in the unit cell, it stores the position of the
        sites.
        """,
    )

    lattice_vectors = Quantity(
        type=np.float64,
        shape=[3, 3],
        unit="angstrom",
        description="""
        Lattice vectors of the simulation cell in cartesian coordinates. The
        last (fastest) index runs over the $x,y,z$ Cartesian coordinates, and the first
        index runs over the 3 lattice vectors.
        """,
    )

    reciprocal_lattice_vectors = Quantity(
        type=np.float64,
        shape=[3, 3],
        unit="1/angstrom",
        description="""
        Reciprocal lattice vectors of the simulation cell, in cartesian coordinates and with the 2 $pi$ pre-factor.
        The first index runs over the $x,y,z$ Cartesian coordinates, and the second index runs
        over the 3 lattice vectors.
        """,
    )

    periodic = Quantity(
        type=bool,
        shape=[3],
        description="""
        Denotes if periodic boundary condition is applied to each of the lattice vectors.'
        """,
    )


class BravaisLattice(MSection):
    """
    Section containing the Bravais lattice metadata.
    """

    m_def = Section(validate=False)

    n_points = Quantity(
        type=np.int32,
        description="""
        Number of Bravais lattice points.
        """,
    )

    points = Quantity(
        type=np.float64,
        shape=["n_points", 3],
        unit="angstrom",
        description="""
        Values of the Bravais lattice points used to obtain the hopping integrals. They are
        sorted from smaller to larger values of the norm.
        """,
    )

    formula_hill = Quantity(
        type=str,
        description="""
        The chemical formula for a structure in Hill form with element symbols followed
        by non-reduced integer chemical proportion numbers. The proportion number is
        omitted if it is 1.
        """,
    )

    system = SubSection(sub_section=System.m_def)


class Model(MSection):
    """
    Section containing the tight-binding model metadata.
    """

    m_def = Section(validate=False)

    bravais_lattice = SubSection(sub_section=BravaisLattice.m_def)

    n_orbitals = Quantity(
        type=np.int32,
        description="""
        Number of projected orbitals.
        """,
    )

    degeneracy_factors = Quantity(
        type=np.int32,
        shape=["n_points"],
        description="""
        Degeneracy of each Bravais lattice point.
        """,
    )

    onsite_energies = Quantity(
        type=np.float64,
        shape=["n_orbitals"],
        unit="eV",
        description="""
        Values of the onsite energies for each orbital.
        """,
    )

    hopping_matrix = Quantity(
        type=np.complex128,
        shape=["n_points", "n_orbitals", "n_orbitals"],
        unit="eV",
        description="""
        Real space hopping matrix for each Bravais lattice point as a matrix of dimension
        (n_orbitals * n_orbitals).

        It is sorted from smaller bravais_lattice.points norm to larger.
        """,
    )
