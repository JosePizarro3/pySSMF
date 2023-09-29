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

import json
import os
import warnings
import numpy as np


class Input:
    _valid_lattice_models = ['linear', 'square', 'honeycomb', 'triangular']

    def __init__(self, **kwargs):
        """Reads the input arguments and stores then in self.data in a JSON file generated
        in the working_directory.
            - 'code': 'SSMF'
            - 'working_directory': working directory in which files are located.
            If 'read_from_input_file' is true:
                - 'input_file': path to the input JSON file in the working directory.
                - 'input_data': read input data from input file.
            Else:
                If 'model_file' is specified:
                    - 'model': path to the tight-binding model file in the working directory.
                    - 'prune_threshold': if 'pruning' is True, it is set to the read or the
                        default 0.01 value.
                If 'lattice_model' is specified:
                    - 'lattice_model': the specific label for the model to be applied (defined
                        in _valid_lattice_models).
                    - 'hoppings': list of values for the hoppings in order from nearest to
                        farthest hopping.
                    - 'n_hoppings': integer specifying how many hoppings are included. Only
                        supported up to 3.
                    - 'n_orbitals': integer specifying how many orbitals are included.

        Then, it stores the input parameters in a JSON in the working directory with the name
        'input_ssmf.json'.
        """
        data = {'code': 'SSMF'}
        if not kwargs.get('working_directory'):
            raise ValueError('Could not find specified the working_directory in the input.')
        data['working_directory'] = kwargs.get('working_directory')
        # Read from input file if provided in the argument
        read_input_file = kwargs.get('read_from_input_file', False)
        if read_input_file:
            data['input_file'] = kwargs.get('input_file')
            input_file = os.path.join(kwargs.get('working_directory'), kwargs.get('input_file'))
            data['input_data'] = self.read_from_file(input_file)
        else:
            # We check whether a model_file has been defined or a model_label
            if kwargs.get('model_file'):
                data['model'] = kwargs.get('model_file')
                if kwargs.get('pruning', False):  # pruning only applies for model_file cases
                    data['prune_threshold'] = kwargs.get('prune_threshold', 0.01)
            elif kwargs.get('lattice_model') in self._valid_lattice_models:
                data['model'] = kwargs.get('lattice_model')
                hoppings = kwargs.get('hoppings', [])
                # We check if 'hoppings' was empty
                if len(hoppings) == 0:
                    n_hoppings = 1
                    n_orbitals = 1
                    hoppings = [[1.0]]
                # Only up to 3 neighbors hoppings supported
                n_hoppings = len(hoppings)
                if n_hoppings > 3:
                    raise ValueError('Maximum n_hoppings models supported is 3. Please, select '
                                     'a smaller number')
                n_orbitals = len(hoppings[0])
                # We check shape for all R point to be (n_orbitals, n_orbitals)
                if not all(np.shape(hop_point) == (n_orbitals, n_orbitals) for hop_point in hoppings):
                    raise ValueError('Dimensions of each hopping matrix do not coincide with '
                                     f'({n_orbitals}, {n_orbitals}).')
                # TODO improve this
                onsite_energies = kwargs.get('onsite_energies', [])
                if len(onsite_energies) == 0:
                    onsite_energies = [0.0] * n_orbitals
                data['hoppings'] = hoppings
                data['onsite_energies'] = onsite_energies
                data['n_hoppings'] = n_hoppings
                data['n_orbitals'] = n_orbitals
            else:
                raise ValueError('Could not find the initial model. Please, check your inputs: '
                                 '1) define `model_file` pointing to your Wannier90 `*_hr.dat` '
                                 'hoppings file, or 2) specify the `lattice_model` to study among the '
                                 f'accepted values {self._valid_lattice_models}.')
        self.data = data
        self.to_json()

    def to_json(self):
        with open(f"{self.data.get('working_directory')}/input_ssmf.json", 'w') as file:
            json.dump(self.data, file, indent=4)

    def read_from_file(self, input_file: str):
        try:
            with open(input_file, 'r') as file:
                input_data = json.load(file)
        except FileNotFoundError:
            raise ValueError(f'Input file {input_file} not found.')
        except json.JSONDecodeError:
            raise ValueError(f'Failed to decode JSON in input file {input_file}.')
        code_name = input_data.get('code', '')
        if code_name != 'SSMF':
            raise ValueError(f'Could not recognize the input JSON file {input_file} as readable by the SSMF code.')
        return input_data
