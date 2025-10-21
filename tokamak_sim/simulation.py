"""Core tokamak simulation implementation."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable, Optional

import numpy as np

from .constants import ALPHA_BIRTH_ENERGY_J, GRAVITY, NEUTRON_ENERGY_J, c
from .fields import FieldSolver
from .particle import Particle
from .species import ALPHA, DEUTERON, ELECTRON, TRITON

DensityProfile = Callable[[float, float], float]
TemperatureProfile = Callable[[float, float], float]


@dataclass
class TokamakSimulation:
    """High level orchestrator for the particle-in-cell tokamak model."""

    major_radius: float = 1.6
    minor_radius: float = 0.45
    toroidal_field: float = 3.5
    dt: float = 1e-8
    max_time: float = 5e-3
    energy_loss_rate: float = 0.5
    fusion_cross_section: float = 1e-21
    max_wall_energy: float = 5e7
    grid_shape: tuple[int, int, int] = (24, 24, 24)
    grid_spacing: tuple[float, float, float] = (0.02, 0.02, 0.02)

    particles: list[Particle] = field(default_factory=list)
    time: float = 0.0
    wall_energy: float = 0.0
    fusion_energy: float = 0.0
    abort_reason: Optional[str] = None

    def __post_init__(self) -> None:
        self.field_solver = FieldSolver(self.grid_shape, self.grid_spacing, self.toroidal_field)

    # ------------------------------------------------------------------
    # Initialisation utilities
    def initialise_plasma(
        self,
        density_profile: DensityProfile,
        temperature_profile: TemperatureProfile,
        macro_particles: int = 5000,
        rng: Optional[np.random.Generator] = None,
    ) -> None:
        rng = np.random.default_rng() if rng is None else rng
        self.particles.clear()
        volume = 2.0 * np.pi * self.major_radius * np.pi * self.minor_radius**2
        max_density = max(density_profile(0.0, 0.0), 1e16)
        while len(self.particles) < macro_particles:
            r = self.minor_radius * rng.random() ** (1.0 / 3.0)
            theta = 2.0 * np.pi * rng.random()
            z = (rng.random() * 2.0 - 1.0) * self.minor_radius
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            density = density_profile(r, z)
            if density <= 0.0:
                continue
            if rng.random() > density / max_density:
                continue
            temperature = temperature_profile(r, z)
            temperature = max(temperature, 1.0)
            ion_species = DEUTERON if rng.random() < 0.5 else TRITON
            vth_ion = np.sqrt(2.0 * temperature * 1.602e-19 / ion_species.mass)
            vth_e = np.sqrt(2.0 * temperature * 1.602e-19 / ELECTRON.mass)
            ion_velocity = rng.normal(0.0, vth_ion / np.sqrt(3.0), size=3)
            elec_velocity = rng.normal(0.0, vth_e / np.sqrt(3.0), size=3)
            pos = np.array([x, y, z])
            self.particles.append(Particle(ion_species, pos, ion_velocity))
            self.particles.append(Particle(ELECTRON, pos.copy(), elec_velocity))
        weight = density_profile(0.0, 0.0) * volume / max(len(self.particles), 1)
        for particle in self.particles:
            particle.weight = weight

    # ------------------------------------------------------------------
    def _deposit_fields(self) -> None:
        self.field_solver.reset()
        for particle in self.particles:
            self.field_solver.deposit(particle.position, particle.species.charge * particle.weight, particle.velocity)
        self.field_solver.solve()

    def _push_particles(self) -> None:
        dt = self.dt
        half = 0.5 * dt
        for particle in self.particles:
            charge = particle.species.charge
            mass = particle.species.mass
            E, B = self.field_solver.gather(particle.position)
            v_minus = particle.velocity + half * (charge / mass) * E
            t = half * (charge / mass) * B
            t_mag2 = float(np.dot(t, t))
            v_prime = v_minus + np.cross(v_minus, t)
            s = 2.0 * t / (1.0 + t_mag2)
            v_plus = v_minus + np.cross(v_prime, s)
            velocity = v_plus + half * (charge / mass) * E
            velocity[2] -= GRAVITY * dt
            particle.velocity = velocity
            particle.position = particle.position + velocity * dt
            particle.update_gamma()

    def _handle_collisions(self) -> None:
        rng = np.random.default_rng()
        angle_std = 0.05
        for particle in self.particles:
            speed = np.linalg.norm(particle.velocity)
            if speed == 0.0:
                continue
            axis = rng.normal(0.0, 1.0, size=3)
            axis = axis - axis.dot(particle.velocity) * particle.velocity / (speed**2)
            axis_norm = np.linalg.norm(axis)
            if axis_norm == 0.0:
                continue
            axis /= axis_norm
            angle = rng.normal(0.0, angle_std)
            particle.velocity = (
                particle.velocity * np.cos(angle)
                + np.cross(axis, particle.velocity) * np.sin(angle)
                + axis * axis.dot(particle.velocity) * (1.0 - np.cos(angle))
            )
            particle.update_gamma()

    def _apply_losses(self) -> None:
        if self.energy_loss_rate <= 0.0:
            return
        factor = np.exp(-self.energy_loss_rate * self.dt / 2.0)
        for particle in self.particles:
            particle.velocity *= factor
            particle.update_gamma()

    def _fusion_events(self) -> None:
        rng = np.random.default_rng()
        volume = np.prod(self.grid_shape) * np.prod(self.grid_spacing)
        rate = self.fusion_cross_section * self.dt / max(volume, 1e-12)
        new_particles: list[Particle] = []
        neutrons = 0.0
        to_remove: set[int] = set()
        for i, pi in enumerate(self.particles):
            if pi.species not in (DEUTERON, TRITON):
                continue
            for j in range(i + 1, len(self.particles)):
                if j in to_remove:
                    continue
                pj = self.particles[j]
                if {pi.species, pj.species} != {DEUTERON, TRITON}:
                    continue
                if rng.random() < rate:
                    direction = rng.normal(0.0, 1.0, size=3)
                    direction /= np.linalg.norm(direction)
                    speed = np.sqrt(2.0 * ALPHA_BIRTH_ENERGY_J / ALPHA.mass)
                    alpha = Particle(ALPHA, 0.5 * (pi.position + pj.position), direction * speed)
                    alpha.weight = 0.5 * (pi.weight + pj.weight)
                    new_particles.append(alpha)
                    neutrons += NEUTRON_ENERGY_J
                    to_remove.add(i)
                    to_remove.add(j)
                    break
        if to_remove:
            self.particles = [p for idx, p in enumerate(self.particles) if idx not in to_remove]
        self.particles.extend(new_particles)
        self.fusion_energy += neutrons

    def _apply_boundary(self) -> None:
        survivors: list[Particle] = []
        for particle in self.particles:
            r = np.linalg.norm(particle.position[:2])
            if r >= self.minor_radius or abs(particle.position[2]) >= self.minor_radius:
                self.wall_energy += particle.kinetic_energy
            else:
                survivors.append(particle)
        self.particles = survivors
        if self.wall_energy > self.max_wall_energy:
            self.abort_reason = "wall energy limit exceeded"

    # ------------------------------------------------------------------
    def step(self) -> None:
        if self.abort_reason is not None:
            return
        self._deposit_fields()
        self._push_particles()
        self._handle_collisions()
        self._apply_losses()
        self._fusion_events()
        self._apply_boundary()
        self.time += self.dt
        if self.time >= self.max_time:
            self.abort_reason = "completed"

    def run(self, controller: Optional["EnergyController"] = None) -> None:
        from .controller import EnergyController

        energy_history: list[float] = []
        time_history: list[float] = []
        while self.abort_reason is None:
            self.step()
            energy = sum(p.kinetic_energy for p in self.particles)
            energy_history.append(energy)
            time_history.append(self.time)
            if controller is not None:
                power = controller.command(energy)
                self._apply_heating(power)
        self.energy_history = np.array(energy_history)
        self.time_history = np.array(time_history)

    def _apply_heating(self, power: float) -> None:
        if power <= 0.0 or not self.particles:
            return
        total_energy = sum(p.kinetic_energy for p in self.particles)
        if total_energy <= 0.0:
            return
        scale = np.sqrt(1.0 + power * self.dt / total_energy)
        for particle in self.particles:
            particle.velocity *= scale
            particle.update_gamma()

    # Convenience accessors
    @property
    def total_energy(self) -> float:
        return sum(p.kinetic_energy for p in self.particles)

    def iter_positions(self) -> Iterable[np.ndarray]:
        for particle in self.particles:
            yield particle.position

    def get_species_indices(self) -> np.ndarray:
        codes = {ELECTRON: 0, DEUTERON: 1, TRITON: 2, ALPHA: 3}
        return np.array([codes.get(p.species, -1) for p in self.particles])
