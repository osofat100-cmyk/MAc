# C++ Port Plan

## Goals
- Reimplement the tokamak particle-in-cell simulation in modern C++ for high performance and future GPU offload.
- Preserve existing physics behaviour so Python <-> C++ results agree within numerical tolerance.
- Provide modular architecture so solvers, controllers, and diagnostics can evolve independently.

## Target Architecture
- `include/sim/` headers and `src/` implementations compiled into a static library plus CLI tools.
- Core modules:
  - `constants.hpp`: physical constants.
  - `species.hpp`: immutable species definitions.
  - `particle.hpp`: structure-of-arrays container for particle state (position, velocity, weight, gamma).
  - `field_solver.hpp`: quasi-static mesh solver (initially CPU Jacobi; later switchable backends via strategy pattern).
  - `controller.hpp`: energy feedback law interface with baseline nonlinear implementation.
  - `simulation.hpp`: orchestrator handling initialise / step / run similar to Python API.
  - `diagnostics.hpp`: streaming energy/wall data and export hooks (CSV/HDF5 planned).
- Math types: Eigen (`Eigen::Vector3d`, `Eigen::Array`) for convenience while we profile SoA kernels.
- Parallelism: CPU OpenMP pragmas around particle pushes and field deposition; future CUDA/Kokkos backends can reuse the SoA layout.
- Build: CMake project targeting C++20, with options for AddressSanitizer and OpenMP toggles.

## Migration Strategy
1. **Parity scaffolding**: set up CMake, basic headers, and stub implementations returning deterministic placeholders. Include GoogleTest harness comparing small runs against Python snapshots.
2. **Particles & energy controller**: port data structures and heating/loss routines; verify kinetic energy calculations using fixtures exported from Python.
3. **Field deposition/solver**: implement Jacobi Poisson solve and gather operations; benchmark vs. Python on identical grid/particle inputs.
4. **Fusion & collisions**: translate stochastic routines, ensuring RNG seeds match Python reference where possible for deterministic tests.
5. **Boundary & diagnostics**: replicate wall loss logic, add structured logging, export to JSON/CSV.
6. **Optimisation passes**: introduce OpenMP, vectorised loops, and optional GPU kernels; profile with VTune/Nsight to guide tuning.

## Validation
- Retain Python harness generating snapshot files (particle states, field arrays, energy history) to compare against C++ outputs.
- Add regression tests for conservation properties (mass/charge, bounded gamma).
- Continuous integration: clang + gcc builds, sanitiser jobs, formatting (clang-format), and unit tests.

### Testing Roadmap
- Add GoogleTest fixtures that load reference snapshots produced by the Python implementation and assert numerical parity for:
  - single-step Boris push (positions, velocities, gamma);
  - field deposition followed by Jacobi solve (E/B components at sample nodes);
  - fusion/collision event statistics over seeded RNG sequences.
- Introduce property tests (rapidcheck or similar) to stress-check invariants such as non-negative weights and bounded gamma.
- Provide optional integration harness that executes a short Python run via subprocess, exports `npz` data, and compares against the compiled binary.

### Immediate Next Steps
1. Port plasma initialisation logic (density/temperature sampling) and expose a builder API for external drivers.
2. Flesh out `FieldSolver` with the full Jacobi solver and deposition/gather routines.
3. Translate energy loss/heating, boundary checks, and controller coupling.
4. Implement fusion and collision sampling using a reproducible RNG (pcg32 family) and add deterministic tests.
5. Wire up CSV/JSON diagnostics plus a simple CLI runner for smoke tests.

## Open Questions
- Select final linear algebra backend (Eigen vs. hand-rolled arrays vs. Kokkos Views) once profiling identifies hotspots.
- Decide on RNG library (PCG, Philox) that maps cleanly to CUDA later.
- Determine data interchange format for Python <-> C++ co-simulations (HDF5 suggested).
