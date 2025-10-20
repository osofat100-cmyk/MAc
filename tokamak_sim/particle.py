"""Particle definitions for the tokamak simulation."""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

from .constants import c
from .species import Species


@dataclass
class Particle:
    """Representation of a macro particle in the simulation."""

    species: Species
    position: np.ndarray
    velocity: np.ndarray
    weight: float = 1.0
    gamma: float = field(init=False)

    def __post_init__(self) -> None:
        self.position = np.asarray(self.position, dtype=float)
        self.velocity = np.asarray(self.velocity, dtype=float)
        self.update_gamma()

    def update_gamma(self) -> None:
        v2 = float(np.dot(self.velocity, self.velocity))
        if v2 >= c ** 2:
            v2 = 0.999999 * c ** 2
        self.gamma = 1.0 / np.sqrt(1.0 - v2 / c ** 2)

    @property
    def kinetic_energy(self) -> float:
        return (self.gamma - 1.0) * self.species.mass * c ** 2 * self.weight

    def momentum(self) -> np.ndarray:
        return self.gamma * self.species.mass * self.velocity * self.weight

    def clone(self) -> "Particle":
        return Particle(
            species=self.species,
            position=self.position.copy(),
            velocity=self.velocity.copy(),
            weight=self.weight,
        )
