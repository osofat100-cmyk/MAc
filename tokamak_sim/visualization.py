"""Matplotlib based 3-D live animation of the tokamak plasma."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from .simulation import TokamakSimulation


@dataclass
class TokamakAnimator:
    """Real-time 3-D scatter animation of the torus and plasma."""

    simulation: TokamakSimulation
    frame_interval: int = 30
    sample_particles: int = 2000
    cmap: str = "viridis"

    def animate(self, controller=None, frames: Optional[int] = None) -> animation.FuncAnimation:
        self.simulation.energy_history = np.array([])
        self.simulation.time_history = np.array([])
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
    lim = 1.2 * self.simulation.minor_radius
    ax.set_xlim((-lim, lim))
    ax.set_ylim((-lim, lim))
    ax.set_zlim((-lim, lim))
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_title("3-D Tokamak Plasma Simulation")
    # larger marker size and slight transparency to improve visibility
    scatter = ax.scatter([], [], [], c=[], cmap=self.cmap, s=12, alpha=0.9, vmin=0, vmax=3)

        torus_u = np.linspace(0, 2 * np.pi, 50)
        torus_v = np.linspace(0, 2 * np.pi, 20)
        u, v = np.meshgrid(torus_u, torus_v)
        R = self.simulation.major_radius
        r = self.simulation.minor_radius
        torus_x = (R + r * np.cos(v)) * np.cos(u) - R
        torus_y = (R + r * np.cos(v)) * np.sin(u)
        torus_z = r * np.sin(v)
        ax.plot_wireframe(torus_x, torus_y, torus_z, color="grey", linewidth=0.3, alpha=0.5)

        # Ensure axes have equal aspect (matplotlib >=3.3)
        try:
            ax.set_box_aspect([1, 1, 1])
        except Exception:
            pass

        def update(frame: int):
            if self.simulation.abort_reason is None:
                self.simulation.step()
                energy = self.simulation.total_energy
                if controller is not None:
                    power = controller.command(energy)
                    self.simulation._apply_heating(power)
            positions = np.array(list(self.simulation.iter_positions()))
            if len(positions) == 0:
                scatter._offsets3d = ([], [], [])
                scatter.set_array(np.array([]))
                return scatter
            idx = np.linspace(0, len(positions) - 1, min(self.sample_particles, len(positions)), dtype=int)
            positions = positions[idx]
            scatter._offsets3d = (positions[:, 0], positions[:, 1], positions[:, 2])
            species = self.simulation.get_species_indices()[idx]
            scatter.set_array(species)
            ax.set_title(f"t = {self.simulation.time*1e3:.2f} ms | Particles = {len(self.simulation.particles)}")
            return scatter

        frame_count = frames or int(self.simulation.max_time / self.simulation.dt)
        ani = animation.FuncAnimation(
            fig,
            update,
            frames=frame_count,
            interval=self.frame_interval,
            blit=False,
            repeat=False,
        )
        plt.tight_layout()
        return ani
