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


class ValidLatticeModels(ABC):
    """Abstract class that defines the valid lattice models covered in this code by specific strings."""

    def __init__(self, logger: logging.Logger = logging.getLogger(__name__)):
        self.logger = logger
        self._valid_lattice_models = [
            "linear",
            "square",
            "honeycomb",
            "triangular",
        ]  # TODO extend this


class Input(ValidLatticeModels):
    def __init__(self, **kwargs):
        """Reads the input arguments and stores then in self.data in a JSON file generated
        in the working_directory.
            - 'logger': the logger where the errors, warnings, etc. will be printed.
            - 'code': 'SSMF'.
            - 'working_directory': working directory where the files are located.
            If 'read_from_input_file' is true:
                - 'input_file': path to the input JSON file in the working directory.
                - 'input_data': read input data from input file.
            Else:
                If 'tb_model_file' is specified:
                    - 'tb_model': path to the tight-binding model file in the working directory.
                    - 'prune_threshold': if 'pruning' is True, it is set to the read or the
                        default 0.01 value.
                If 'lattice_model' is specified:
                    - 'tb_model': the specific label for the model to be applied (defined
                        in self._valid_lattice_models).
                    - 'hoppings': list of values for the hoppings in order from nearest to
                        farthest hopping.
                    - 'n_hoppings': integer specifying how many hoppings are included. Only
                        supported up to 3.
                    - 'n_orbitals': integer specifying how many orbitals are included.

        Then, it stores the input parameters in a JSON in the working directory with the name
        'input_ssmf.json'.
        """
        super().__init__()
        # Initializing data
        data = {"code": "SSMF"}

        # Check working_directory and stores it in data
        if not kwargs.get("working_directory"):
            self.logger.error(
                "Could not find specified the working_directory in the input."
            )
            return
        data["working_directory"] = kwargs.get("working_directory")

        # Read from input file if provided in the argument
        read_input_file = kwargs.get("read_from_input_file", False)
        if read_input_file:
            data["input_file"] = kwargs.get("input_file")
            input_file = os.path.join(
                kwargs.get("working_directory"), kwargs.get("input_file")
            )
            data["input_data"] = self.read_from_file(input_file)
        else:
            # We check whether a model_file has been defined or a model_label
            if kwargs.get("tb_model_file"):
                data["tb_model"] = kwargs.get("tb_model_file")
                # pruning only applies for tb_model_file cases
                if kwargs.get("pruning", False):
                    data["prune_threshold"] = kwargs.get("prune_threshold", 0.01)
            elif kwargs.get("lattice_model") in self._valid_lattice_models:
                data["tb_model"] = kwargs.get("lattice_model")
                hoppings = kwargs.get("hoppings", [])
                # We check if 'hoppings' was empty
                if len(hoppings) == 0:
                    n_hoppings = 1
                    n_orbitals = 1
                    hoppings = [[1.0]]
                # Only up to 3 neighbors hoppings supported
                n_hoppings = len(hoppings)
                if n_hoppings > 3:
                    self.logger.error(
                        "Maximum n_hoppings models supported is 3. Please, select "
                        "a smaller number."
                    )
                    return
                n_orbitals = len(hoppings[0])
                # We check shape for all R point to be (n_orbitals, n_orbitals)
                if not all(
                    np.shape(hop_point) == (n_orbitals, n_orbitals)
                    for hop_point in hoppings
                ):
                    self.logger.error(
                        "Dimensions of each hopping matrix do not coincide with"
                        "(n_orbitals, n_orbitals).",
                        data={"n_orbitals": n_orbitals},
                    )
                # TODO improve this
                onsite_energies = kwargs.get("onsite_energies", [])
                if len(onsite_energies) == 0:
                    onsite_energies = [0.0] * n_orbitals
                data["hoppings"] = hoppings
                data["onsite_energies"] = onsite_energies
                data["n_hoppings"] = n_hoppings
                data["n_orbitals"] = n_orbitals
            else:
                self.logger.error(
                    "Could not find the initial model. Please, check your inputs: "
                    "1) define `model_file` pointing to your Wannier90 `*_hr.dat` "
                    "hoppings file, or 2) specify the `lattice_model` to study among the "
                    "accepted values.",
                    data={"lattice_model": self._valid_lattice_models},
                )

            # KGrids
            # For band structure calculations
            data["n_k_path"] = kwargs.get("n_k_path", 90)
            # For full_bz diagonalization
            data["k_grid"] = kwargs.get("k_grid", [1, 1, 1])

            # Plotting arguments
            data["plot_hoppings"] = kwargs.get("plot_hoppings", False)
            data["plot_bands"] = kwargs.get("plot_bands", False)
        self.data = data
        self.to_json()

    def to_json(self):
        with open(f"{self.data.get('working_directory')}/input_ssmf.json", "w") as file:
            json.dump(self.data, file, indent=4)

    def read_from_file(self, input_file: str):
        try:
            with open(input_file, "r") as file:
                input_data = json.load(file)
        except FileNotFoundError:
            self.logger.error(
                "Input file not found.",
                extra={"input_file": input_file},
            )
        except json.JSONDecodeError:
            self.logger.error(
                "Failed to decode JSON in input file.",
                extra={"input_file": input_file},
            )
        code_name = input_data.get("code", "")
        if code_name != "SSMF":
            self.logger.error(
                "Could not recognize the input JSON file as readable by the SSMF code.",
                extra={"input_file": input_file},
            )
        return input_data
