"""Core physics primitives for the tokamak simulation.

This module implements lightweight data structures for particle species,
macro-particles, and a finite-difference field solver operating on a 3D grid.
The implementation focuses on clarity and extensibility rather than raw
performance so that the classes can serve as educational references or the
foundation for higher fidelity solvers.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

# Physical constants (SI units)
QE = 1.602176634e-19  # Elementary charge (C)
ME = 9.10938356e-31   # Electron mass (kg)
MP = 1.67262192369e-27  # Proton mass (kg)
C = 2.99792458e8      # Speed of light (m/s)
MU0 = 4e-7 * np.pi    # Vacuum permeability (H/m)
EPS0 = 1.0 / (MU0 * C ** 2)  # Vacuum permittivity (F/m)
G = 9.81              # Gravitational acceleration (m/s^2)


@dataclass(frozen=True)
class Species:
    """A charged particle species.

    Attributes
    ----------
    name:
        Descriptive name of the species.
    mass:
        Particle mass in kilograms.
    charge:
        Particle charge in coulombs.
    """

    name: str
    mass: float
    charge: float


@dataclass
class Particle:
    """A macro-particle used by the particle-in-cell solver."""

    species: Species
    position: np.ndarray
    velocity: np.ndarray

    def __post_init__(self) -> None:
        self.position = np.asarray(self.position, dtype=float)
        self.velocity = np.asarray(self.velocity, dtype=float)
        if self.position.shape != (3,) or self.velocity.shape != (3,):
            raise ValueError("position and velocity must be 3-vectors")
        self.gamma = self._compute_gamma()

    def kinetic_energy(self) -> float:
        """Return the kinetic energy in joules."""

        return (self.gamma - 1.0) * self.species.mass * C ** 2

    def momentum(self) -> np.ndarray:
        """Return the relativistic momentum vector."""

        return self.gamma * self.species.mass * self.velocity

    def advance(self, displacement: np.ndarray) -> None:
        """Translate the particle in-place."""

        self.position += displacement

    def set_velocity(self, new_velocity: np.ndarray) -> None:
        """Assign a new velocity and update the relativistic factor."""

        self.velocity = np.asarray(new_velocity, dtype=float)
        self.gamma = self._compute_gamma()

    def _compute_gamma(self) -> float:
        vmag2 = float(np.dot(self.velocity, self.velocity))
        if vmag2 >= C ** 2:
            # Clamp to avoid numerical blow-up during pathological steps
            return 1e6
        return 1.0 / np.sqrt(1.0 - vmag2 / C ** 2)


class FieldSolver:
    """Electrostatic and magnetostatic solver on a regular 3D grid."""

    def __init__(self, grid_shape: Iterable[int], grid_spacing: Iterable[float], *,
                 toroidal_field: float) -> None:
        self.grid_shape = tuple(int(v) for v in grid_shape)
        self.grid_spacing = tuple(float(v) for v in grid_spacing)
        if len(self.grid_shape) != 3 or len(self.grid_spacing) != 3:
            raise ValueError("grid_shape and grid_spacing must be length 3")
        self.toroidal_field = float(toroidal_field)

        self.charge_density = np.zeros(self.grid_shape, dtype=float)
        self.current_density = np.zeros((*self.grid_shape, 3), dtype=float)
        self.electric_field = np.zeros((*self.grid_shape, 3), dtype=float)
        self.magnetic_field = np.zeros((*self.grid_shape, 3), dtype=float)

    # ------------------------------------------------------------------
    # Charge and current deposition
    # ------------------------------------------------------------------
    def reset(self) -> None:
        self.charge_density.fill(0.0)
        self.current_density.fill(0.0)

    def deposit(self, particle: Particle) -> None:
        """Deposit charge and current to the nearest grid cell."""

        ix, iy, iz = self._cell_index(particle.position)
        if ix is None:
            return
        self.charge_density[ix, iy, iz] += particle.species.charge
        self.current_density[ix, iy, iz] += particle.species.charge * particle.velocity

    # ------------------------------------------------------------------
    # Field solution
    # ------------------------------------------------------------------
    def solve(self, iterations: int = 40) -> None:
        """Compute electric and magnetic fields from deposited sources."""

        potential = np.zeros(self.grid_shape, dtype=float)
        dx, dy, dz = self.grid_spacing

        for _ in range(iterations):
            potential[1:-1, 1:-1, 1:-1] = (
                potential[:-2, 1:-1, 1:-1]
                + potential[2:, 1:-1, 1:-1]
                + potential[1:-1, :-2, 1:-1]
                + potential[1:-1, 2:, 1:-1]
                + potential[1:-1, 1:-1, :-2]
                + potential[1:-1, 1:-1, 2:]
                - self.charge_density[1:-1, 1:-1, 1:-1] * dx ** 2 / EPS0
            ) / 6.0

        ex = (np.roll(potential, -1, axis=0) - np.roll(potential, 1, axis=0)) / (2 * dx)
        ey = (np.roll(potential, -1, axis=1) - np.roll(potential, 1, axis=1)) / (2 * dy)
        ez = (np.roll(potential, -1, axis=2) - np.roll(potential, 1, axis=2)) / (2 * dz)
        self.electric_field[..., 0] = -ex
        self.electric_field[..., 1] = -ey
        self.electric_field[..., 2] = -ez

        curl_x = (
            (np.roll(self.current_density[..., 2], -1, axis=1) -
             np.roll(self.current_density[..., 2], 1, axis=1)) / (2 * dy)
            - (np.roll(self.current_density[..., 1], -1, axis=2) -
               np.roll(self.current_density[..., 1], 1, axis=2)) / (2 * dz)
        )
        curl_y = (
            (np.roll(self.current_density[..., 0], -1, axis=2) -
             np.roll(self.current_density[..., 0], 1, axis=2)) / (2 * dz)
            - (np.roll(self.current_density[..., 2], -1, axis=0) -
               np.roll(self.current_density[..., 2], 1, axis=0)) / (2 * dx)
        )
        curl_z = (
            (np.roll(self.current_density[..., 1], -1, axis=0) -
             np.roll(self.current_density[..., 1], 1, axis=0)) / (2 * dx)
            - (np.roll(self.current_density[..., 0], -1, axis=1) -
               np.roll(self.current_density[..., 0], 1, axis=1)) / (2 * dy)
        )

        self.magnetic_field[..., 0] = MU0 * curl_x * dx
        self.magnetic_field[..., 1] = MU0 * curl_y * dy + self.toroidal_field
        self.magnetic_field[..., 2] = MU0 * curl_z * dz

    # ------------------------------------------------------------------
    # Field gathering
    # ------------------------------------------------------------------
    def gather(self, particle: Particle) -> tuple[np.ndarray, np.ndarray]:
        ix, iy, iz = self._cell_index(particle.position)
        if ix is None:
            return np.zeros(3), np.zeros(3)
        return self.electric_field[ix, iy, iz], self.magnetic_field[ix, iy, iz]

    # ------------------------------------------------------------------
    def _cell_index(self, position: np.ndarray) -> tuple[int | None, int, int] | tuple[None, None, None]:
        nx, ny, nz = self.grid_shape
        dx, dy, dz = self.grid_spacing
        ix = int(np.floor(position[0] / dx + nx / 2))
        iy = int(np.floor(position[1] / dy + ny / 2))
        iz = int(np.floor(position[2] / dz + nz / 2))
        if 0 <= ix < nx and 0 <= iy < ny and 0 <= iz < nz:
            return ix, iy, iz
        return None, None, None


# Convenience factory for common species -------------------------------------

def default_species() -> dict[str, Species]:
    """Return a dictionary with commonly used species."""

    return {
        "electron": Species("electron", ME, -QE),
        "deuteron": Species("deuteron", 2 * MP, QE),
        "triton": Species("triton", 3 * MP, QE),
        "alpha": Species("alpha", 4 * MP, 2 * QE),
    }
