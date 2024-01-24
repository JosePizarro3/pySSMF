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

matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def plot_hopping_matrices(
    matrices: np.ndarray, max_value: float = 1.0, display_count: int = 4
):
    """
    Plots the matrices with a slider to scroll through number of Bravais points ($N_R$).

    Args:
        matrices (np.array): The matrices to be plotted.
        max_value (float, optional): The maximum value for colormap normalization. Defaults to 1.0.
        display_count (int, optional): Number of matrices to be displayed at once. Defaults to 4.
    """

    total_matrices = matrices.shape[0]

    fig, axes = plt.subplots(
        nrows=1, ncols=display_count, figsize=(display_count * 3, 3)
    )
    if display_count == 1:
        axes = [axes]

    # Create a function to update the displayed matrices
    def update(val):
        idx = int(slider.val)
        for i, ax in enumerate(axes):
            ax.clear()
            if idx + i < total_matrices:
                ax.imshow(matrices[idx + i], cmap="inferno", vmin=0, vmax=max_value)
                ax.set_title(f"$N_R$ = {idx + i}")
            ax.axis("off")

        fig.canvas.draw_idle()

    # Create the slider
    ax_slider = plt.axes([0.2, 0.01, 0.65, 0.03], facecolor="lightgoldenrodyellow")
    slider = Slider(
        ax_slider,
        "Start $N_R$",
        0,
        total_matrices - display_count,
        valinit=0,
        valstep=1,
    )
    slider.on_changed(update)

    # Initial plot
    update(0)

    # Add a colorbar to the figure, adjust the position for better spacing
    cbar_ax = fig.add_axes([0.93, 0.15, 0.02, 0.7])
    fig.colorbar(
        plt.cm.ScalarMappable(
            cmap="inferno", norm=plt.Normalize(vmin=0, vmax=max_value)
        ),
        cax=cbar_ax,
    )

    plt.tight_layout(
        rect=[0, 0.05, 0.9, 1]
    )  # Adjust the right bound to 0.9 to provide space for the colorbar
    plt.show()


def plot_band_structure(eigenvalues, tb_hamiltonian, special_points=None):
    """
    Plots the band structure of a Hamiltonian.

    Args:
        eigenvalues (np.ndarray): Eigenvalues of the Hamiltonian matrix of shape (Nk, Norb).
        special_points (dict, optional): Dictionary of special points and their labels.
    """

    num_bands = eigenvalues.shape[1]

    # Create a figure
    plt.figure(figsize=(8, 6))

    # Iterate over eigenvalues for each band
    for band_idx in range(num_bands):
        plt.plot(
            np.arange(len(eigenvalues)),
            eigenvalues[:, band_idx],
            label=f"Band {band_idx + 1}",
        )

    # Customize x-axis labeling based on special points
    if special_points is not None:
        x_labels = [""] * (len(special_points))  # Initialize labels with empty strings
        x_ticks = []

        i = 0
        for key, val in special_points.items():
            if np.where(np.all(val == tb_hamiltonian.k_path.kpts, axis=1))[0].size > 0:
                index = np.where(np.all(val == tb_hamiltonian.k_path.kpts, axis=1))[0][
                    0
                ]
                x_labels[i] = key
                x_ticks.append(index)
                i += 1
        plt.xticks(x_ticks, x_labels)

    plt.xlim(0, len(eigenvalues) - 1)
    plt.xlabel("k-points")
    plt.ylabel("Energy (eV)")
    plt.title("Band Structure")
    plt.legend()
    plt.grid(True)

    plt.show()


def plot_dos(energies, orbital_dos, total_dos):
    """
    Plots the density of states (DOS) of a tight-binding Hamiltonian.

    Args:
        energies: the energies at which the DOS is evaluated.
        orbital_dos: the orbital-resolved DOS.
        total_dos: the total DOS.
    """
    # Create a figure
    plt.figure(figsize=(8, 6))
    plt.plot(energies, total_dos, label="Total DOS", color="k", linewidth=3.5)

    # Plot orbital-resolved DOS
    for i, orb_dos in enumerate(orbital_dos):
        plt.plot(energies, orb_dos, label=f"Orbital {i + 1}")

    plt.xlabel("Energy (eV)")
    plt.ylabel("Density of States (DOS)")
    plt.legend()
    plt.show()
