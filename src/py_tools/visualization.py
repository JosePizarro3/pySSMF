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
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def plot_hopping_matrices(matrices: np.ndarray, max_value: float = 1.0, display_count: int = 4):
    """
    Plot the matrices with a slider to scroll through number of Bravais points ($N_R$).

    Args:
        matrices (np.array): The matrices to be plotted.
        max_value (float, optional): The maximum value for colormap normalization. Defaults to 1.0.
        display_count (int, optional): Number of matrices to be displayed at once. Defaults to 4.
    """

    total_matrices = matrices.shape[0]

    fig, axes = plt.subplots(nrows=1, ncols=display_count, figsize=(display_count*3, 3))
    if display_count == 1:
        axes = [axes]

    # Create a function to update the displayed matrices
    def update(val):
        idx = int(slider.val)
        for i, ax in enumerate(axes):
            ax.clear()
            if idx + i < total_matrices:
                ax.imshow(matrices[idx + i], cmap='inferno', vmin=0, vmax=max_value)
                ax.set_title(f'$N_R$ = {idx + i}')
            ax.axis('off')

        fig.canvas.draw_idle()

    # Create the slider
    ax_slider = plt.axes([0.2, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'Start $N_R$', 0, total_matrices - display_count, valinit=0, valstep=1)
    slider.on_changed(update)

    # Initial plot
    update(0)

    # Add a colorbar to the figure, adjust the position for better spacing
    cbar_ax = fig.add_axes([0.93, 0.15, 0.02, 0.7])
    fig.colorbar(plt.cm.ScalarMappable(cmap='inferno', norm=plt.Normalize(vmin=0, vmax=max_value)), cax=cbar_ax)

    plt.tight_layout(rect=[0, 0.05, 0.9, 1])  # Adjust the right bound to 0.9 to provide space for the colorbar
    plt.show()
