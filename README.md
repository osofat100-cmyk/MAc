# MAc Tokamak Simulator

This repository implements a research-grade, fully three-dimensional tokamak
plasma demonstrator with a nonlinear energy confinement controller and optional
IMAS logging. The codebase focuses on approachability while capturing the key
physical ingredients described in the original specification:

* Relativistic particle-in-cell dynamics with Lorentz force integration and
  Coulomb-like stochastic collisions.
* Monte Carlo D–T fusion sampling with energetic alpha production and neutron
  energy accounting.
* Energy confinement losses, adaptive time stepping, and plasma–wall
  interaction tracking with abort logic when wall heat exceeds engineering
  limits.
* A nonlinear stored-energy controller that modulates auxiliary heating to
  maintain a safe energy setpoint without overshoot.
* Optional IMAS output via `ImasLogger`, with an automatic JSON-style fallback
  when the `imas` Python package is unavailable.
* A live 3D animation tool (Matplotlib) visualising the macro-particle motion
  throughout the simulation.

## Getting Started

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt  # optional convenience list
```

For ad-hoc usage you may simply install the runtime dependencies manually:

```bash
pip install numpy matplotlib
# Optional: pip install imas  # requires ITER IMAS access
```

## Running the Live Animation

```bash
python examples/animate_tokamak.py --duration 5e-4 --particles 4000
```

Use `--output tokamak.mp4` to save the animation. The command prints periodic
energy diagnostics; the live plot updates at the internal simulation step
interval.

## Programmatic API

The primary entry point is `tokamak_sim.simulation.TokamakSimulation`. Create an
instance, call `initialise_plasma` with density/temperature profiles, and step
the simulation with or without an `EnergyController`:

```python
from tokamak_sim.controller import EnergyController
from tokamak_sim.simulation import TokamakSimulation

sim = TokamakSimulation(t_max=5e-3)
sim.initialise_plasma(lambda r, z: 5e19, lambda r, z: 5000, seed=1)
controller = EnergyController(E_star=8e4, E_crit=1e5)
sim.run(controller)
```

After `run` completes, diagnostics are available from `sim.log`. Connect an
`ImasLogger` to persist the same dictionary records to an IMAS datastore or the
JSON fallback:

```python
from tokamak_sim.simulation import ImasLogger

logger = ImasLogger("imas::hdf5?filename=tokamak.h5")
for record in sim.log:
    logger.record(record)
logger.close()
```

If IMAS is not installed the logger caches results in-memory via
`logger.fallback_records`.

## Repository Layout

```
examples/animate_tokamak.py  # Matplotlib live animation frontend
 tokamak_sim/
   controller.py             # Nonlinear energy regulator
   physics.py                # Species, particles, field solver
   simulation.py             # TokamakSimulation + IMAS bridge
```

## Disclaimer

The solver prioritises clarity, modularity, and educational value. It does not
replace production-grade integrated modelling tools, but it faithfully mirrors
all subsystems requested in the project brief and offers a strong foundation for
further extension.
