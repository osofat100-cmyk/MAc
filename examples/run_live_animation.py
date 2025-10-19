"""Entry point to launch the 3-D tokamak animation."""

from __future__ import annotations

import numpy as np

from tokamak_sim import EnergyController, TokamakAnimator, TokamakSimulation


def main() -> None:
    simulation = TokamakSimulation(max_time=5e-3, dt=2e-7, grid_shape=(20, 20, 20))

    def density_profile(r: float, z: float) -> float:
        return 5e19 * max(0.0, 1.0 - (r**2 + z**2) / simulation.minor_radius**2)

    def temperature_profile(r: float, z: float) -> float:
        return 5e3 * max(0.1, 1.0 - 0.5 * (r**2 + z**2) / simulation.minor_radius**2)

    simulation.initialise_plasma(density_profile, temperature_profile, macro_particles=4000)
    controller = EnergyController(target_energy=8e4, critical_energy=1e5)

    animator = TokamakAnimator(simulation, frame_interval=50, sample_particles=1000)
    animation = animator.animate(controller=controller, frames=400)
    animation.save("tokamak.mp4", writer="ffmpeg")


if __name__ == "__main__":
    main()
