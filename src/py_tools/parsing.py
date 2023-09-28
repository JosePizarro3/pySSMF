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
import logging
import json

# NOMAD functionalities for parsing the tight-binding calculation files
from nomad.parsing.file_parser import TextParser, Quantity
from nomad.units import ureg
# NOMAD schema
from . import get_files, System, BravaisLattice, Model


re_n = r'[\n\r]'


class WOutParser(TextParser):
    """Parses the text from the `*.wout` Wannier90 file. It finds (using regex) the values
    for the atomic structure.
    """
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        structure_quantities = [
            Quantity('labels', r'\|\s*([A-Z][a-z]*)', repeats=True),
            Quantity('positions', r'\|\s*([\-\d\.]+)\s*([\-\d\.]+)\s*([\-\d\.]+)', repeats=True)
        ]

        self._quantities = [
            Quantity(
                'lattice_vectors', r'\s*a_\d\s*([\d\-\s\.]+)', repeats=True),
            Quantity(
                'structure', rf'(\s*Fractional Coordinate[\s\S]+?)(?:{re_n}\s*(PROJECTIONS|K-POINT GRID))',
                repeats=False, sub_parser=TextParser(quantities=structure_quantities)),
        ]


class HrParser(TextParser):
    """Parses the `*_hr.dat` Wannier90 file. It finds (using regex) the degeneracy factors
    at the begining of the file and the list of hoppings.
    """
    def __init__(self):
        super().__init__(None)

    def init_quantities(self):
        self._quantities = [
            Quantity('degeneracy_factors', r'\s*written on[\s\w]*:\d*:\d*\s*([\d\s]+)'),
            Quantity('hoppings', r'\s*([-\d\s.]+)', repeats=False)
        ]


class MinimalWannier90Parser():
    def __init__(self):
        self.wout_parser = WOutParser()
        self.hr_parser = HrParser()

    def init_parser(self):
        """Initializes the parser attributes.
        """
        self.wout_parser.mainfile = self.mainfile
        self.wout_parser.logger = self.logger
        self.hr_parser.mainfile = self.filepath
        self.hr_parser.logger = self.logger

    def parse_system(self, model: Model):
        """Parses the system metadata.

        Args:
            model (Model): section Model to store the system information.
        """
        sec_system = model.m_create(BravaisLattice).m_create(System)

        structure = self.wout_parser.get('structure')
        if structure is None:
            self.logger.error('Error parsing the structure from .wout')
            return
        if self.wout_parser.get('lattice_vectors', []):
            lattice_vectors = np.vstack(self.wout_parser.get('lattice_vectors', [])[-3:])
            sec_system.lattice_vectors = lattice_vectors * ureg.angstrom
            sec_system.periodic = [True, True, True]
        sec_system.labels = structure.get('labels')
        if structure.get('positions') is not None:
            sec_system.positions = structure.get('positions') * ureg.angstrom

    def parse_hoppings(self, model: Model):
        """Parses the hoppings metadata.

        Args:
            model (Model): section Model to store the hoppings information.
        """
        bravais_lattice = model.bravais_lattice

        deg_factors = self.hr_parser.get('degeneracy_factors', [])
        full_hoppings = self.hr_parser.get('hoppings', [])
        if deg_factors is not None and full_hoppings is not None:
            n_orbitals = deg_factors[0]
            n_points = deg_factors[1]
            wann90_hops = np.reshape(full_hoppings, (n_points, n_orbitals, n_orbitals, 7))
            # Parse BravaisLattice
            bravais_lattice.n_points = n_points
            bravais_lattice_points = wann90_hops[:, 0, 0, :3]
            bravais_lattice_points = bravais_lattice_points @ bravais_lattice.system.lattice_vectors.magnitude
            sorted_indices = np.argsort(np.linalg.norm(bravais_lattice_points, axis=1))
            bravais_lattice_points_sorted = bravais_lattice_points[sorted_indices]
            bravais_lattice.points = bravais_lattice_points_sorted
            # Parse degeneracy factors
            deg_factors = deg_factors[2:]
            deg_factors_sorted = deg_factors[sorted_indices]
            model.degeneracy_factors = deg_factors_sorted
            model.n_orbitals = n_orbitals
            hops = wann90_hops[:, :, :, 5] + 1j * wann90_hops[:, :, :, 6]
            hops_sorted = hops[sorted_indices]
            onsite_energies = []
            for morb in range(n_orbitals):
                onsite_energies.append(hops_sorted[0][morb][morb].real)
                hops_sorted[0][morb][morb] = 0.0 + hops_sorted[0][morb][morb].imag
            model.onsite_energies = onsite_energies * ureg.eV
            model.hopping_matrix = hops_sorted * ureg.eV

    def parse(self, filepath: str, model: Model, logger: logging.Logger = None):
        """Parses the system, Bravais lattice and hoppings information and stores it in the
        section Model.

        Args:
            filepath (str): path to the file `*_hr.dat` to be parsed.
            model (Model): section Model to store metadata.
        """
        basename = os.path.basename(filepath)  # Getting filepath for *_hr.dat file
        wout_files = get_files('*.wout', filepath, basename)
        if len(wout_files) > 1:
            logger.warning('Multiple `*.wout` files found; we will parse the last one.')
        mainfile = wout_files[-1]  # Path to *.wout file

        self.filepath = filepath
        self.mainfile = mainfile
        self.logger = logging.getLogger(__name__) if logger is None else logger

        self.init_parser()

        self.parse_system(model)
        self.parse_hoppings(model)


class MinimalTBStudioParser():
    def __init__(self):
        pass

    def init_parser(self):
        pass

    def parse(self, filepath: str, model: Model, logger: logging.Logger = None):  # pylint: disable=unused-argument
        # TODO implement this parsing with the help of @MohamadNakhaee
        self.filepath = filepath
        self.logger = logging.getLogger(__name__) if logger is None else logger
        self.init_parser()


def read_input(input_filepath: str, logger: logging.Logger = None) -> dict:
    """Reads the SSMF code-specific JSON input file and returns a dictionary with the read
    parameters.

    Args:
        input_filepath (str): input JSON file.

    Returns:
        dict: dictionary of input parameters for running SSMF.
    """
    logger = logging.getLogger(__name__) if logger is None else logger
    try:
        input_data = {}
        with open(input_filepath, 'r') as file:
            input_data = json.load(file)
    except FileNotFoundError:
        logger.error(f'Input file {input_filepath} not found.')
    except json.JSONDecodeError:
        logger.error(f'Failed to decode JSON in input file {input_filepath}.')

    code_name = input_data.get('code', '')
    if code_name == 'SSMF':
        return input_data
    else:
        logger.error(f'Could not recognize the input JSON file {input_filepath} as readable by the SSMF code.')
