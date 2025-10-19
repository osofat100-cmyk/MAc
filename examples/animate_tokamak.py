"""Run a short 3D tokamak simulation and render a live animation."""
from __future__ import annotations

import argparse
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from tokamak_sim.controller import EnergyController
from tokamak_sim.simulation import TokamakSimulation


def constant_profile(value: float):
    def profile(r: float, z: float) -> float:
        return value
    return profile


def create_simulation(args: argparse.Namespace) -> tuple[TokamakSimulation, EnergyController]:
    sim = TokamakSimulation(t_max=args.duration)
    density = constant_profile(args.density)
    temperature = constant_profile(args.temperature)
    sim.initialise_plasma(density, temperature, particle_count=args.particles, seed=args.seed)
    controller = EnergyController(E_star=args.target_energy, E_crit=args.critical_energy,
                                  Kp=args.kp, alpha=args.alpha, max_power=args.max_power)
    return sim, controller


def animate(args: argparse.Namespace) -> None:
    sim, controller = create_simulation(args)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_xlim(-sim.minor_radius, sim.minor_radius)
    ax.set_ylim(-sim.minor_radius, sim.minor_radius)
    ax.set_zlim(-sim.minor_radius, sim.minor_radius)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_zlabel("z (m)")
    ax.set_title("Tokamak Plasma Macro-Particle Motion")

    scatter = ax.scatter([], [], [], s=5, alpha=0.6)

    def update(_frame: int):
        sim._step(controller)
        if not sim.particles:
            return scatter,
        positions = np.array([p.position for p in sim.particles])
        scatter._offsets3d = (positions[:, 0], positions[:, 1], positions[:, 2])
        ax.set_title(f"t = {sim.time * 1e3:.2f} ms | Particles = {len(sim.particles)}")
        return scatter,

    frames = int(args.duration / sim.dt)
    ani = animation.FuncAnimation(fig, update, frames=frames, interval=1, blit=False)
    if args.output:
        ani.save(args.output, fps=30)
    else:
        plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration", type=float, default=5e-4, help="Simulation duration (s)")
    parser.add_argument("--density", type=float, default=5e19, help="Core density (m^-3)")
    parser.add_argument("--temperature", type=float, default=5000, help="Core temperature (eV)")
    parser.add_argument("--particles", type=int, default=3000, help="Number of macro-particles")
    parser.add_argument("--seed", type=int, default=1, help="RNG seed for reproducibility")
    parser.add_argument("--target-energy", type=float, default=8e4, help="Energy setpoint (J)")
    parser.add_argument("--critical-energy", type=float, default=1e5, help="Critical energy limit (J)")
    parser.add_argument("--kp", type=float, default=0.2, help="Proportional gain")
    parser.add_argument("--alpha", type=float, default=0.8, help="Nonlinear gain")
    parser.add_argument("--max-power", type=float, default=5e7, help="Maximum heating power (W)")
    parser.add_argument("--output", type=str, help="Write animation to a file instead of showing a window")
    return parser.parse_args()


if __name__ == "__main__":
    animate(parse_args())
