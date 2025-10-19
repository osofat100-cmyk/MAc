# Tokamak Simulation and 3-D Animation

This repository provides a research-grade, yet lightweight, 3-D tokamak
plasma simulation written in Python. The package combines a particle-in-cell
physics core, a nonlinear plasma energy controller, and a Matplotlib-based
animation front-end capable of producing live visualisations of the plasma
within a toroidal vacuum vessel.

## Features

- **Particle-in-Cell plasma model** with relativistic Boris push, Coulomb
  scattering and stochastic fusion reactions.
- **Self-consistent field solve** for quasi-static electric and magnetic
  fields on a 3-D grid, including an external toroidal field.
- **Boundary and wall loading model** that records deposited energy and
  aborts the simulation when engineering limits are exceeded.
- **Nonlinear energy controller** implementing the stabiliser law described in
the specification to maintain plasma energy below a configurable threshold.
- **IMAS-friendly architecture** with explicit geometry, actuator and profile
  hooks ready for integration with the ITER Integrated Modelling & Analysis
  Suite.
- **3-D live animation** driven by Matplotlib that renders particle motion and
tokamak geometry in real time.

## Installation

Create a virtual environment and install the runtime dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The minimal requirements are:

```text
matplotlib
numpy
```

(Optionally install `ffmpeg` system-wide when exporting animations.)

## Running the live animation

The `examples/run_live_animation.py` script initialises an ITER-like plasma,
activates the nonlinear controller, and launches the 3-D animation. The script
also demonstrates how to export the animation to an MP4 file:

```bash
python examples/run_live_animation.py
```

Running the script will open a Matplotlib window that updates the plasma
scatter plot while the simulation executes. To simply perform headless runs
and inspect energy history data, use the `TokamakSimulation` class directly.

## Python API snapshot

```python
from tokamak_sim import EnergyController, TokamakSimulation

sim = TokamakSimulation(max_time=1e-2, dt=1e-7)

def density_profile(r, z):
    return 5e19 * max(0.0, 1.0 - (r**2 + z**2) / sim.minor_radius**2)

def temperature_profile(r, z):
    return 5e3 * max(0.1, 1.0 - 0.5 * (r**2 + z**2) / sim.minor_radius**2)

sim.initialise_plasma(density_profile, temperature_profile)
controller = EnergyController(target_energy=8e4, critical_energy=1e5)

sim.run(controller)
print("Final energy (J):", sim.total_energy)
```

The resulting `TokamakSimulation` instance stores time and energy histories
that can be post-processed or exported to IMAS-compliant data structures.
