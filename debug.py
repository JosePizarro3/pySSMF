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

# General imports
import logging
import os
import numpy as np

from src.py_tools import Model, TBHamiltonian, Pruner, MinimalWannier90Parser, MinimalTBStudioParser, read_input
# Plotting
from src.py_tools.visualization import plot_hopping_matrices, plot_band_structure


# Set up logger
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# Read the input file data
working_dir = os.path.join(os.path.dirname(__file__), 'work_dir')
input_data = read_input(f'{working_dir}/input.json', logger)

# Creating Model for using parsing with NOMAD functionalities in schema.py and parsing.py
model = Model()
# And calling the specific Parser
file = input_data.get('file_tbmodel', '')
basename = os.path.basename(file)
if '_hr.dat' in basename:
    MinimalWannier90Parser().parse(file, model, logger)
elif '.tbm' in basename:
    MinimalTBStudioParser().parse(file, model, logger)


# Pruning hoppings
prune_threshold = input_data.get('prune_threshold', 0.01)
pruner = Pruner(model)
pruner.prune_by_threshold(prune_threshold, logger)
# plot_hopping_matrices(pruner.hopping_matrix_norms / pruner.max_value, 1.0)

# Diagonalizing for band structures
tb_hamiltonian = TBHamiltonian(model)
# tb_hamiltonian.write_to_hdf5(f'{working_dir}/SSMF_output.h5')
special_points = tb_hamiltonian.k_path.special_points
_, eigenvalues = tb_hamiltonian.diagonalize(tb_hamiltonian.kpoints)
plot_band_structure(eigenvalues, tb_hamiltonian, special_points)
