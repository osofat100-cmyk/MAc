from __future__ import annotations

"""Matplotlib-based visualisation helpers for `TokamakSimulation`."""

from dataclasses import dataclass
from typing import Any, Optional, cast

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from .simulation import TokamakSimulation


@dataclass
class TokamakAnimator:
    """Render a 3-D scatter animation of simulation particle positions.

    The `steps_per_frame` parameter lets the animation advance multiple
    simulation steps between redraws, useful when the simulation uses
    a very small `dt` for stability.
    """

    simulation: TokamakSimulation
    frame_interval: int = 30
    sample_particles: int = 2000
    cmap: str = "viridis"
    steps_per_frame: int = 1

    def animate(
        self,
        controller: Optional[Any] = None,
        frames: Optional[int] = None,
    ) -> animation.FuncAnimation:
        """Create and return a matplotlib `FuncAnimation` for the simulation."""

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        ax_any = cast(Any, ax)

        lim = 1.2 * float(getattr(self.simulation, "minor_radius", 1.0))
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_title("3-D Tokamak Plasma Simulation")

        scatter = ax_any.scatter([], [], [], c=[], cmap=self.cmap, s=8, alpha=0.9)

        # Optional torus wireframe for reference
        torus_u = np.linspace(0, 2 * np.pi, 50)
        torus_v = np.linspace(0, 2 * np.pi, 20)
        u, v = np.meshgrid(torus_u, torus_v)
        major_radius = float(getattr(self.simulation, "major_radius", 1.0))
        minor_radius = float(getattr(self.simulation, "minor_radius", 0.3))
        torus_x = (major_radius + minor_radius * np.cos(v)) * np.cos(u) - major_radius
        torus_y = (major_radius + minor_radius * np.cos(v)) * np.sin(u)
        torus_z = minor_radius * np.sin(v)
        ax.plot_wireframe(torus_x, torus_y, torus_z, color="grey", linewidth=0.3, alpha=0.5)

        try:
            ax.set_box_aspect([1, 1, 1])
        except Exception:
            # Older matplotlib versions do not provide set_box_aspect
            pass

        def controller_power(energy: float) -> float:
            if controller is None:
                return 0.0

            for candidate in ("command", "compute_command"):
                method_obj = getattr(controller, candidate, None)
                if callable(method_obj):
                    try:
                        result = method_obj(energy)
                        return float(cast(Any, result))
                    except Exception:
                        return 0.0
            return 0.0

        heating = getattr(self.simulation, "_apply_heating", None)

        def update(_frame: int):
            nonlocal scatter

            steps_this_frame = max(1, int(self.steps_per_frame))
            for _ in range(steps_this_frame):
                if getattr(self.simulation, "abort_reason", None) is not None:
                    break
                self.simulation.step()
                total_energy = float(getattr(self.simulation, "total_energy", 0.0))
                power = controller_power(total_energy)
                if callable(heating):
                    try:
                        heating(power)
                    except Exception:
                        pass

            positions = np.asarray(list(self.simulation.iter_positions()), dtype=float)
            if positions.size == 0:
                try:
                    scatter._offsets3d = ([], [], [])
                except Exception:
                    scatter.remove()
                    scatter = ax_any.scatter([], [], [], c=[], cmap=self.cmap, s=8, alpha=0.9)
                scatter.set_array(np.array([]))
                ax.set_title("3-D Tokamak Plasma Simulation (no particles)")
                return (scatter,)

            nplot = min(self.sample_particles, len(positions))
            if nplot < len(positions):
                idx = np.linspace(0, len(positions) - 1, nplot, dtype=int)
                plot_pos = positions[idx]
            else:
                plot_pos = positions

            try:
                scatter._offsets3d = (plot_pos[:, 0], plot_pos[:, 1], plot_pos[:, 2])
            except Exception:
                # Some matplotlib combos require recreating the collection
                scatter.remove()
                scatter = ax_any.scatter(
                    plot_pos[:, 0],
                    plot_pos[:, 1],
                    plot_pos[:, 2],
                    cmap=self.cmap,
                    s=8,
                    alpha=0.9,
                )

            try:
                species_all = np.asarray(self.simulation.get_species_indices(), dtype=float)
                species = species_all[: len(plot_pos)]
                if species.size < len(plot_pos):
                    species = np.pad(species, (0, len(plot_pos) - species.size))
                scatter.set_array(species)
            except Exception:
                scatter.set_array(np.zeros(len(plot_pos), dtype=float))

            particle_count = len(getattr(self.simulation, "particles", []))
            current_time = float(getattr(self.simulation, "time", 0.0))
            ax.set_title(
                f"t = {current_time * 1e3:.2f} ms | Particles = {particle_count}"
            )
            return (scatter,)

        if frames is None:
            dt = float(getattr(self.simulation, "dt", 1e-3)) or 1e-3
            max_time = float(getattr(self.simulation, "max_time", dt))
            computed = int(max_time / max(dt, 1e-12))
            frame_count = max(1, computed)
        else:
            frame_count = max(1, int(frames))

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
