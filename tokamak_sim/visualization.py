from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Any, cast

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from .simulation import TokamakSimulation


@dataclass
class TokamakAnimator:
    simulation: TokamakSimulation
    frame_interval: int = 30
    sample_particles: int = 2000
    cmap: str = "viridis"

    def animate(self, controller: Optional[Any] = None, frames: Optional[int] = None) -> animation.FuncAnimation:
        """Minimal 3D animation for TokamakSimulation."""
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        ax_any = cast(Any, ax)

        lim = 1.2 * getattr(self.simulation, "minor_radius", 1.0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)

        def controller_power(energy: float) -> float:
            if controller is None:
                return 0.0
            if hasattr(controller, "command") and callable(getattr(controller, "command")):
                return float(controller.command(energy))
            if hasattr(controller, "compute_command") and callable(getattr(controller, "compute_command")):
                return float(controller.compute_command(energy))
            return 0.0

        def update(_frame: int):
            if getattr(self.simulation, "abort_reason", None) is None:
                self.simulation.step()
                energy = float(getattr(self.simulation, "total_energy", 0.0))
                power = controller_power(energy)
                try:
                    self.simulation._apply_heating(power)
                except Exception:
                    pass

            positions = np.array(list(self.simulation.iter_positions()))
            if positions.size == 0:
                return (ax_any.scatter([], [], [], c=[], cmap=self.cmap, s=8),)

            nplot = min(self.sample_particles, len(positions))
            if nplot < len(positions):
                idx = np.linspace(0, len(positions) - 1, nplot, dtype=int)
                plot_pos = positions[idx]
            else:
                plot_pos = positions

            try:
                species_all = np.array(self.simulation.get_species_indices())
                species = species_all[: len(plot_pos)]
            except Exception:
                species = np.zeros(len(plot_pos))

            scat = ax_any.scatter(plot_pos[:, 0], plot_pos[:, 1], plot_pos[:, 2], c=species, cmap=self.cmap, s=8)
            return (scat,)

        frame_count = frames or int(getattr(self.simulation, "max_time", 1.0) / max(getattr(self.simulation, "dt", 1e-3), 1e-12))
        ani = animation.FuncAnimation(fig, update, frames=frame_count, interval=self.frame_interval, blit=False, repeat=False)
        plt.tight_layout()
        return ani
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Any, cast

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from .simulation import TokamakSimulation


@dataclass
class TokamakAnimator:
    simulation: TokamakSimulation
    frame_interval: int = 30
    sample_particles: int = 2000
    cmap: str = "viridis"

    def animate(self, controller: Optional[Any] = None, frames: Optional[int] = None) -> animation.FuncAnimation:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        ax_any = cast(Any, ax)

        lim = 1.2 * getattr(self.simulation, "minor_radius", 1.0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)

        def controller_power(energy: float) -> float:
            if controller is None:
                return 0.0
            if hasattr(controller, "command") and callable(getattr(controller, "command")):
                return float(controller.command(energy))
            if hasattr(controller, "compute_command") and callable(getattr(controller, "compute_command")):
                return float(controller.compute_command(energy))
            return 0.0

        def update(_frame: int):
            if getattr(self.simulation, "abort_reason", None) is None:
                self.simulation.step()
                energy = float(getattr(self.simulation, "total_energy", 0.0))
                power = controller_power(energy)
                try:
                    self.simulation._apply_heating(power)
                except Exception:
                    pass

            positions = np.array(list(self.simulation.iter_positions()))
            if positions.size == 0:
                scat = ax_any.scatter([], [], [], c=[], cmap=self.cmap, s=8)
                return (scat,)

            nplot = min(self.sample_particles, len(positions))
            if nplot < len(positions):
                idx = np.linspace(0, len(positions) - 1, nplot, dtype=int)
                plot_pos = positions[idx]
            else:
                plot_pos = positions

            try:
                species_all = np.array(self.simulation.get_species_indices())
                species = species_all[: len(plot_pos)]
            except Exception:
                species = np.zeros(len(plot_pos))

            scat = ax_any.scatter(plot_pos[:, 0], plot_pos[:, 1], plot_pos[:, 2], c=species, cmap=self.cmap, s=8)
            return (scat,)

        frame_count = frames or int(getattr(self.simulation, "max_time", 1.0) / max(getattr(self.simulation, "dt", 1e-3), 1e-12))
        ani = animation.FuncAnimation(fig, update, frames=frame_count, interval=self.frame_interval, blit=False, repeat=False)
        plt.tight_layout()
        return ani
"""Minimal Matplotlib 3D animator for TokamakSimulation.

This file is kept intentionally compact and type-checker friendly.  It
uses casts when calling the 3-D scatter API to silence Pylance/Pyright
argument signature checks and keeps a single `TokamakAnimator` class.
"""

from __future__ import annotations

from dataclasses import dataclass
"""Matplotlib 3-D animator for TokamakSimulation.

Single, compact implementation to display particle positions in 3-D.
This file is intentionally conservative with type-checker friendly
casts when calling the 3-D scatter API.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Any, cast

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from .simulation import TokamakSimulation


@dataclass
class TokamakAnimator:
    simulation: TokamakSimulation
    frame_interval: int = 30
    sample_particles: int = 2000
    cmap: str = "viridis"

    def animate(self, controller: Optional[Any] = None, frames: Optional[int] = None) -> animation.FuncAnimation:
        """Create a matplotlib FuncAnimation for the simulation."""

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        ax_any = cast(Any, ax)

        lim = 1.2 * getattr(self.simulation, "minor_radius", 1.0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_title("3-D Tokamak Plasma Simulation")

        # Torus wireframe
        torus_u = np.linspace(0, 2 * np.pi, 50)
        torus_v = np.linspace(0, 2 * np.pi, 20)
        u, v = np.meshgrid(torus_u, torus_v)
        R = getattr(self.simulation, "major_radius", 1.0)
        r = getattr(self.simulation, "minor_radius", 0.3)
        torus_x = (R + r * np.cos(v)) * np.cos(u) - R
        torus_y = (R + r * np.cos(v)) * np.sin(u)
        torus_z = r * np.sin(v)
        ax.plot_wireframe(torus_x, torus_y, torus_z, color="grey", linewidth=0.3, alpha=0.5)

        def controller_power(energy: float) -> float:
            if controller is None:
                return 0.0
            if hasattr(controller, "command") and callable(getattr(controller, "command")):
                return float(controller.command(energy))
            if hasattr(controller, "compute_command") and callable(getattr(controller, "compute_command")):
                return float(controller.compute_command(energy))
            return 0.0

        def update(_frame: int):
            # Step simulation and optionally apply controller
            if getattr(self.simulation, "abort_reason", None) is None:
                self.simulation.step()
                energy = float(getattr(self.simulation, "total_energy", 0.0))
                power = controller_power(energy)
                try:
                    self.simulation._apply_heating(power)
                except Exception:
                    pass

            positions = np.array(list(self.simulation.iter_positions()))

            if positions.size == 0:
                scat = ax_any.scatter([], [], [], c=[], cmap=self.cmap, s=8, alpha=0.9)
                return (scat,)

            nplot = min(self.sample_particles, len(positions))
            if nplot < len(positions):
                idx = np.linspace(0, len(positions) - 1, nplot, dtype=int)
                plot_pos = positions[idx]
            else:
                plot_pos = positions

            try:
                species_all = np.array(self.simulation.get_species_indices())
                species = species_all[: len(plot_pos)]
            except Exception:
                species = np.zeros(len(plot_pos))

            scat = ax_any.scatter(plot_pos[:, 0], plot_pos[:, 1], plot_pos[:, 2], c=species, cmap=self.cmap, s=8, alpha=0.9)
            ax.set_title(f"t = {getattr(self.simulation, 'time', 0.0)*1e3:.2f} ms | Particles = {len(self.simulation.particles)}")
            return (scat,)

        frame_count = frames or int(getattr(self.simulation, "max_time", 1.0) / max(getattr(self.simulation, "dt", 1e-3), 1e-12))
        ani = animation.FuncAnimation(fig, update, frames=frame_count, interval=self.frame_interval, blit=False, repeat=False)
        plt.tight_layout()
        return ani


import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from .simulation import TokamakSimulation


@dataclass
class TokamakAnimator:
    simulation: TokamakSimulation
    frame_interval: int = 30
    sample_particles: int = 2000
    cmap: str = "viridis"

    def animate(self, controller: Optional[object] = None, frames: Optional[int] = None) -> animation.FuncAnimation:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")

        lim = 1.2 * getattr(self.simulation, "minor_radius", 1.0)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_title("3-D Tokamak Plasma Simulation")

        # Torus wireframe for reference
        torus_u = np.linspace(0, 2 * np.pi, 50)
        torus_v = np.linspace(0, 2 * np.pi, 20)
        u, v = np.meshgrid(torus_u, torus_v)
        R = getattr(self.simulation, "major_radius", 1.0)
        r = getattr(self.simulation, "minor_radius", 0.3)
        torus_x = (R + r * np.cos(v)) * np.cos(u) - R
        torus_y = (R + r * np.cos(v)) * np.sin(u)
        torus_z = r * np.sin(v)
        ax.plot_wireframe(torus_x, torus_y, torus_z, color="grey", linewidth=0.3, alpha=0.5)

        def controller_power(energy: float) -> float:
            if controller is None:
                return 0.0
            if hasattr(controller, "command") and callable(getattr(controller, "command")):
                return float(controller.command(energy))
            if hasattr(controller, "compute_command") and callable(getattr(controller, "compute_command")):
                return float(controller.compute_command(energy))
            return 0.0

        def update(_frame: int):
            # Advance simulation and optionally apply controller action
            if getattr(self.simulation, "abort_reason", None) is None:
                self.simulation.step()
                energy = float(getattr(self.simulation, "total_energy", 0.0))
                power = controller_power(energy)
                try:
                    self.simulation._apply_heating(power)
                except Exception:
                    pass

            positions = np.array(list(self.simulation.iter_positions()))

            # If no particles, return empty collection
            if positions.size == 0:
                scat = ax.scatter([], [], [], c=[], cmap=self.cmap, s=8, alpha=0.9)
                return (scat,)

            # Sample for plotting performance
            nplot = min(self.sample_particles, len(positions))
            if nplot < len(positions):
                idx = np.linspace(0, len(positions) - 1, nplot, dtype=int)
                plot_pos = positions[idx]
            else:
                plot_pos = positions

            # Color by species if available
            try:
                species_all = np.array(self.simulation.get_species_indices())
                species = species_all[: len(plot_pos)]
            except Exception:
                species = np.zeros(len(plot_pos))

            scat = ax.scatter(plot_pos[:, 0], plot_pos[:, 1], plot_pos[:, 2], c=species, cmap=self.cmap, s=8, alpha=0.9)
            ax.set_title(f"t = {getattr(self.simulation, 'time', 0.0)*1e3:.2f} ms | Particles = {len(self.simulation.particles)}")
            return (scat,)

        frame_count = frames or int(getattr(self.simulation, "max_time", 1.0) / max(getattr(self.simulation, "dt", 1e-3), 1e-12))
        ani = animation.FuncAnimation(fig, update, frames=frame_count, interval=self.frame_interval, blit=False, repeat=False)
        plt.tight_layout()
        return ani
"""Matplotlib based 3-D live animation of the tokamak plasma."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Any, cast

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
        # Prepare figure and axes
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")
        lim = 1.2 * self.simulation.minor_radius
        """Matplotlib based 3-D live animation of the tokamak plasma.

        This module provides a single clean TokamakAnimator class that runs a
        matplotlib 3-D scatter animation of particles produced by
        TokamakSimulation. It is intentionally small and defensive: it will
        work with or without an external controller object and tolerates an
        empty particle set.
        """

        from __future__ import annotations

        from dataclasses import dataclass
        from typing import Optional, Any, cast

        import matplotlib.pyplot as plt
        from matplotlib import animation
        import numpy as np

        from .simulation import TokamakSimulation


        @dataclass
        class TokamakAnimator:
            """Real-time 3-D scatter animation of the torus and plasma.

            Attributes:
                simulation: TokamakSimulation instance to animate.
                frame_interval: milliseconds between frames.
                sample_particles: maximum number of particles to plot each frame.
                cmap: matplotlib colormap name to color by species index.
            """

            simulation: TokamakSimulation
            frame_interval: int = 30
            sample_particles: int = 2000
            cmap: str = "viridis"

            def animate(self, controller: Optional[object] = None, frames: Optional[int] = None) -> animation.FuncAnimation:
                """Create and return a matplotlib FuncAnimation for the simulation.

                The optional `controller` may provide either a `command(energy)`
                or `compute_command(energy)` method; if present it will be used
                to compute heating power applied each step.
                """

                # Prepare figure and axes
                fig = plt.figure(figsize=(8, 6))
                ax = fig.add_subplot(111, projection="3d")
                lim = 1.2 * getattr(self.simulation, "minor_radius", 1.0)
                ax.set_xlim(-lim, lim)
                ax.set_ylim(-lim, lim)
                ax.set_zlim(-lim, lim)
                ax.set_xlabel("x [m]")
                ax.set_ylabel("y [m]")
                ax.set_zlabel("z [m]")
                ax.set_title("3-D Tokamak Plasma Simulation")

            # Empty 3D scatter initialization
            scatter = cast(Any, ax).scatter(np.array([]), np.array([]), np.array([]), c=np.array([]), cmap=self.cmap, s=8, alpha=0.9)

                # Torus wireframe for reference (shifted so centre is near origin)
                torus_u = np.linspace(0, 2 * np.pi, 50)
                torus_v = np.linspace(0, 2 * np.pi, 20)
                u, v = np.meshgrid(torus_u, torus_v)
                R = getattr(self.simulation, "major_radius", 1.0)
                r = getattr(self.simulation, "minor_radius", 0.3)
                torus_x = (R + r * np.cos(v)) * np.cos(u) - R
                torus_y = (R + r * np.cos(v)) * np.sin(u)
                torus_z = r * np.sin(v)
                ax.plot_wireframe(torus_x, torus_y, torus_z, color="grey", linewidth=0.3, alpha=0.5)

                # Try to set equal aspect ratio where supported
                try:
                    ax.set_box_aspect([1, 1, 1])
                except Exception:
                    pass

                def controller_power(energy: float) -> float:
                    if controller is None:
                        return 0.0
                    if hasattr(controller, "command") and callable(getattr(controller, "command")):
                        return float(controller.command(energy))
                    if hasattr(controller, "compute_command") and callable(getattr(controller, "compute_command")):
                        return float(controller.compute_command(energy))
                    return 0.0

                def update(_frame: int):
                    nonlocal scatter

                    # Advance simulation and optionally apply controller heating
                    if getattr(self.simulation, "abort_reason", None) is None:
                        self.simulation.step()
                        energy = float(getattr(self.simulation, "total_energy", 0.0))
                        power = controller_power(energy)
                        try:
                            # use internal heating hook if available
                            self.simulation._apply_heating(power)
                        except Exception:
                            pass

                    positions = np.array(list(self.simulation.iter_positions()))
                    if positions.size == 0:
                        # recreate an empty scatter
                        try:
                            scatter.remove()
                        except Exception:
                            pass
                        scatter = cast(Any, ax).scatter(np.array([]), np.array([]), np.array([]), c=np.array([]), cmap=self.cmap, s=8, alpha=0.9)
                        return (scatter,)

                    # choose an index subset to plot for performance
                    nplot = min(self.sample_particles, len(positions))
                    if nplot < len(positions):
                        idx = np.linspace(0, len(positions) - 1, nplot, dtype=int)
                        plot_pos = positions[idx]
                    else:
                        plot_pos = positions

                    # color by species index if available
                    try:
                        species_all = np.array(self.simulation.get_species_indices())
                        species = species_all[: len(plot_pos)]
                    except Exception:
                        species = np.zeros(len(plot_pos))

                    # recreate scatter with new points/colors (robust across matplotlib versions)
                try:
                    scatter.remove()
                except Exception:
                    pass
                scatter = cast(Any, ax).scatter(plot_pos[:, 0], plot_pos[:, 1], plot_pos[:, 2], c=species, cmap=self.cmap, s=8, alpha=0.9)

                    ax.set_title(f"t = {getattr(self.simulation, 'time', 0.0)*1e3:.2f} ms | Particles = {len(self.simulation.particles)}")
                    return (scatter,)

                frame_count = frames or int(getattr(self.simulation, "max_time", 1.0) / max(getattr(self.simulation, "dt", 1e-3), 1e-12))
                ani = animation.FuncAnimation(fig, update, frames=frame_count, interval=self.frame_interval, blit=False, repeat=False)
                plt.tight_layout()
                return ani
        ax.set_zlabel("z [m]")
        ax.set_title("3-D Tokamak Plasma Simulation")

        # Empty 3D scatter initialization
        scatter = ax.scatter([], [], [], c=[], cmap=self.cmap, s=8, alpha=0.9)

        # Torus wireframe for reference (shifted so centre is near origin)
        torus_u = np.linspace(0, 2 * np.pi, 50)
        torus_v = np.linspace(0, 2 * np.pi, 20)
        u, v = np.meshgrid(torus_u, torus_v)
        R = getattr(self.simulation, "major_radius", 1.0)
        r = getattr(self.simulation, "minor_radius", 0.3)
        torus_x = (R + r * np.cos(v)) * np.cos(u) - R
        torus_y = (R + r * np.cos(v)) * np.sin(u)
        torus_z = r * np.sin(v)
        ax.plot_wireframe(torus_x, torus_y, torus_z, color="grey", linewidth=0.3, alpha=0.5)

        # Try to set equal aspect ratio where supported
        try:
            ax.set_box_aspect([1, 1, 1])
        except Exception:
            pass

        def controller_power(energy: float) -> float:
            if controller is None:
                return 0.0
            if hasattr(controller, "command") and callable(getattr(controller, "command")):
                return float(controller.command(energy))
            if hasattr(controller, "compute_command") and callable(getattr(controller, "compute_command")):
                return float(controller.compute_command(energy))
            return 0.0

        def update(_frame: int):
            # Advance simulation and optionally apply controller heating
            if getattr(self.simulation, "abort_reason", None) is None:
                self.simulation.step()
                energy = float(getattr(self.simulation, "total_energy", 0.0))
                power = controller_power(energy)
                try:
                    # use internal heating hook if available
                    self.simulation._apply_heating(power)
                except Exception:
                    pass

            positions = np.array(list(self.simulation.iter_positions()))
            if positions.size == 0:
                # clear scatter
                try:
                    scatter._offsets3d = ([], [], [])
                except Exception:
                    pass
                scatter.set_array(np.array([]))
                return (scatter,)

            # choose an index subset to plot for performance
            nplot = min(self.sample_particles, len(positions))
            if nplot < len(positions):
                idx = np.linspace(0, len(positions) - 1, nplot, dtype=int)
                plot_pos = positions[idx]
            else:
                plot_pos = positions

            # update 3-D scatter offsets
            try:
                scatter._offsets3d = (plot_pos[:, 0], plot_pos[:, 1], plot_pos[:, 2])
            except Exception:
                # fallback: re-create the scatter collection (slower)
                scatter.remove()
                c = ax.scatter(plot_pos[:, 0], plot_pos[:, 1], plot_pos[:, 2], c=np.zeros(len(plot_pos)), cmap=self.cmap, s=8, alpha=0.9)
                return (c,)

            # color by species index if available
            try:
                species_all = np.array(self.simulation.get_species_indices())
                species = species_all[: len(plot_pos)]
                scatter.set_array(species)
            except Exception:
                scatter.set_array(np.zeros(len(plot_pos)))

            ax.set_title(f"t = {getattr(self.simulation, 'time', 0.0)*1e3:.2f} ms | Particles = {len(self.simulation.particles)}")
            return (scatter,)

        frame_count = frames or int(getattr(self.simulation, "max_time", 1.0) / max(getattr(self.simulation, "dt", 1e-3), 1e-12))
        ani = animation.FuncAnimation(fig, update, frames=frame_count, interval=self.frame_interval, blit=False, repeat=False)
        plt.tight_layout()
        return ani
