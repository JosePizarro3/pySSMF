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

import os
from glob import glob
import h5py
from typing import Union
import numpy as np


def get_files(pattern: str, filepath: str, stripname: str = "", deep: bool = True):
    """Get files following the `pattern` with respect to the file `stripname` (usually this
    being the mainfile of the given parser) up to / down from the `filepath` (`deep=True` going
    down, `deep=False` up)

    Args:
        pattern (str): targeted pattern to be found
        filepath (str): filepath to start the search
        stripname (str, optional): name with respect to which do the search. Defaults to ''.
        deep (bool, optional): boolean setting the path in the folders to scan (up or down). Defaults to True.

    Returns:
        list: List of found files.
    """
    for _ in range(10):
        filenames = glob(f"{os.path.dirname(filepath)}/{pattern}")
        pattern = os.path.join("**" if deep else "..", pattern)
        if filenames:
            break

    if len(filenames) > 1:
        # filter files that match
        suffix = os.path.basename(filepath).strip(stripname)
        matches = [f for f in filenames if suffix in f]
        filenames = matches if matches else filenames

    filenames = [f for f in filenames if os.access(f, os.F_OK)]
    return filenames


def extract_hdf5_dataset(
    data: Union[h5py.Dataset, h5py.Group],
    default: Union[bool, int, float, np.ndarray, None] = None,
):
    """Extracts the hdf5 dataset and returns it as its Python native type. It can also return
    a default value if specified and data is not a `h5py.Dataset` object.

    Args:
        data (Union[h5py.Dataset, h5py.Group]): Input h5py data.
        default (Union[bool, int, float, np.ndarray, None], optional): Output default value.
            Defaults to None.
    """
    if isinstance(data, h5py.Dataset):
        return data[()]
    else:
        return default
