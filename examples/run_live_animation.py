"""Entry point to launch the 3-D tokamak animation."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from tokamak_sim import EnergyController, TokamakAnimator, TokamakSimulation


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--frames",
        type=int,
        default=400,
        help="Number of animation frames to render (default: %(default)s)",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save the animation to an MP4 file instead of displaying it live.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tokamak.mp4"),
        help="Output path used when --save is specified (default: %(default)s)",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    simulation = TokamakSimulation(max_time=5e-3, dt=2e-7, grid_shape=(20, 20, 20))

    def density_profile(r: float, z: float) -> float:
        return 5e19 * max(0.0, 1.0 - (r**2 + z**2) / simulation.minor_radius**2)

    def temperature_profile(r: float, z: float) -> float:
        return 5e3 * max(0.1, 1.0 - 0.5 * (r**2 + z**2) / simulation.minor_radius**2)

    simulation.initialise_plasma(density_profile, temperature_profile, macro_particles=4000)
    controller = EnergyController(target_energy=8e4, critical_energy=1e5)

    animator = TokamakAnimator(simulation, frame_interval=50, sample_particles=1000)
    animation = animator.animate(controller=controller, frames=args.frames)

    if args.save:
        animation.save(args.output, writer="ffmpeg")
    else:
        plt.show()


if __name__ == "__main__":
    main()
