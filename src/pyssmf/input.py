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

from abc import ABC
import logging
import json
import os
import numpy as np

from . import LOGGER


class Input:
    def __init__(self, **kwargs):
        """
        Reads the input arguments and stores then in a dictionary called `data` and in a
        JSON file, `input_ssmf.json` generated in the working_directory.

        Attributes:
            data (dict): A dictionary that stores input arguments. The keys and their
            corresponding values are as follows:
                - 'code' (str): Always set to 'pySSMF'.
                - 'working_directory' (str): The directory where the files are located.
                - 'input_file' (str): Path to the input JSON file in the working directory.
                Only present if 'read_from_input_file' is True.
                - 'lattice_model' (str): the specific lattice model to be calculated. See
                `valid_lattice_models` for accepted values.
                - 'hoppings' (list): List of values for the hoppings in order from nearest to
                farthest hopping. Only present if 'lattice_model' is specified.
                - 'n_hoppings' (int): Integer specifying how many hoppings are included.
                Only supported up to 3. Only present if 'lattice_model' is specified.
                - 'n_orbitals' (int): Integer specifying how many orbitals are included.
                Only present if 'lattice_model' is specified.
        """
        super().__init__()
        # List of covered lattice toy models
        # TODO extend this list
        _valid_lattices = [
            'linear',
        ]

        # Reading input from an `input.json` file
        read_input_file = kwargs.get('read_from_input_file', False)
        if read_input_file:
            input_file = os.path.join(
                kwargs.get('working_directory'), kwargs.get('input_file')
            )
            data = self.read_from_file(input_file)
            data['input_file'] = kwargs.get('input_file')
            self.data = data
            self.to_json()
            return

        # If `input_file` is not specified, we populate `data` with the passed arguments instead
        # Initializing `data`
        data = {'code': 'pySSMF'}

        # Check working_directory and stores it in data
        if not kwargs.get('working_directory'):
            raise KeyError(
                'Could not find specified the working_directory in the input.'
            )
        data['working_directory'] = kwargs.get('working_directory')

        # Lattice models details
        lattice_model_id = kwargs.get('lattice_model', '')
        if lattice_model_id not in _valid_lattices:
            raise ValueError(f'{lattice_model_id} is not a valid lattice model.')
        data['lattice_model'] = lattice_model_id
        # We check if 'hoppings' was empty
        hoppings = kwargs.get('hoppings', [])
        n_hoppings = len(hoppings)
        if n_hoppings == 0:
            n_hoppings = 1
            n_orbitals = 1
            hoppings = [[1.0]]
            LOGGER.warning(
                'Argument `hoppings` was empty, so we consider a nearest neighbor, single-orbital model'
            )
        # Only up to 3 neighbors hoppings supported
        n_hoppings = len(hoppings)
        if n_hoppings > 3:
            raise ValueError(
                'Maximum n_hoppings models supported is 3. Please, select '
                'a smaller number.'
            )
        # We check shape for all Wigner-Seitz points hoppings to be (n_orbitals, n_orbitals)
        n_orbitals = len(hoppings[0])
        if not all(
            np.shape(hop_point) == (n_orbitals, n_orbitals) for hop_point in hoppings
        ):
            raise ValueError(
                'Dimensions of each hopping matrix do not coincide with'
                '(n_orbitals, n_orbitals).',
                data={'n_orbitals': n_orbitals},
            )
        # Extracting `onsite_energies`, `hoppings`
        onsite_energies = kwargs.get('onsite_energies', [])
        if len(onsite_energies) == 0:
            onsite_energies = [0.0] * n_orbitals
            LOGGER.warning(
                'Attribute `onsite_energies` was empty, so we consider all zeros with the dimensions of `(n_orbitals)`.'
            )
        data['hoppings'] = hoppings
        data['onsite_energies'] = onsite_energies
        data['n_hoppings'] = n_hoppings
        data['n_orbitals'] = n_orbitals

        # KGrids
        # For band structure calculations
        data['n_k_path'] = kwargs.get('n_k_path', 90)
        # For full_bz diagonalization
        data['k_grid'] = kwargs.get('k_grid', [1, 1, 1])

        # Plotting arguments
        data['plot_hoppings'] = kwargs.get('plot_hoppings', False)
        data['plot_bands'] = kwargs.get('plot_bands', False)
        # DOS calculation and plotting
        data['dos'] = kwargs.get('dos', False)
        data['dos_gaussian_width'] = kwargs.get('dos_gaussian_width', 0.1)
        data['dos_delta_energy'] = kwargs.get('dos_delta_energy', 0.01)
        # Nominal number of electrons
        data['n_electrons'] = kwargs.get('n_electrons', 1)
        self.data = data
        self.to_json()

    def to_json(self) -> None:
        """
        Stores the input data in a JSON file in the working directory.
        """
        with open(f"{self.data.get('working_directory')}/input_ssmf.json", 'w') as file:
            json.dump(self.data, file, indent=4)

    def read_from_file(self, input_file: str) -> dict:
        """
        Reads the input data from a JSON file provided it is a pySSMF code input file.

        Args:
            input_file (str): path to the input JSON file in the working directory.

        Returns:
            (dict): dictionary with the input data read from the JSON input file.
        """
        try:
            with open(input_file, 'r') as file:
                input_data = json.load(file)
        except (FileNotFoundError, json.JSONDecodeError):
            raise FileNotFoundError(
                'Input file not found or failed to decode JSON input file.',
                extra={'input_file': input_file},
            )
        code_name = input_data.get('code', '')
        if code_name != 'pySSMF':
            raise ValueError(
                'Could not recognize the input JSON file as readable by the pySSMF code.',
                extra={'input_file': input_file},
            )
        return input_data
