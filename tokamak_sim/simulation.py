"""High-level 3D tokamak simulation with optional IMAS logging."""
from __future__ import annotations

from dataclasses import dataclass, field
import importlib.util
import math
from typing import Callable, List, Sequence

import numpy as np

from .controller import EnergyController
from .physics import (C, G, Particle, default_species,
                      FieldSolver)

# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def _has_imas() -> bool:
    return importlib.util.find_spec("imas") is not None


# ---------------------------------------------------------------------------
# Tokamak geometry helpers
# ---------------------------------------------------------------------------


def cylindrical_radius(position: np.ndarray, major_radius: float) -> float:
    x, y, _ = position
    return math.sqrt((x + major_radius) ** 2 + y ** 2) - major_radius


@dataclass
class BoundaryResult:
    lost_particles: int = 0
    deposited_energy: float = 0.0
    divertor_energy: float = 0.0
    wall_energy: float = 0.0


@dataclass
class TokamakSimulation:
    """A lightweight particle-in-cell tokamak model.

    Parameters
    ----------
    major_radius:
        Tokamak major radius (m).
    minor_radius:
        Tokamak minor radius (m).
    toroidal_field:
        Nominal toroidal magnetic field at the magnetic axis (T).
    grid_shape/grid_spacing:
        Resolution of the field solver grid.
    dt:
        Initial time step (s).
    energy_loss_rate:
        Exponential energy loss coefficient (1/s) representing transport and
        radiation.
    fusion_reactivity:
        Simplified reactivity used for Monte Carlo fusion sampling (m^3/s).
    max_wall_energy:
        Abort threshold for energy deposited on plasma-facing components (J).
    """

    major_radius: float = 1.6
    minor_radius: float = 0.45
    toroidal_field: float = 3.5
    grid_shape: Sequence[int] = (32, 32, 32)
    grid_spacing: Sequence[float] = (0.02, 0.02, 0.02)
    dt: float = 1e-8
    t_max: float = 1e-2
    energy_loss_rate: float = 0.5
    fusion_reactivity: float = 1e-21
    max_wall_energy: float = 5e7

    particles: List[Particle] = field(default_factory=list)
    time: float = 0.0

    def __post_init__(self) -> None:
        self.field_solver = FieldSolver(self.grid_shape, self.grid_spacing,
                                        toroidal_field=self.toroidal_field)
        self.species = default_species()
        self.total_neutron_energy = 0.0
        self.cumulative_wall_energy = 0.0
        self.abort_triggered = False
        self.log: list[dict[str, float]] = []

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------
    def initialise_plasma(
        self,
        density_profile: Callable[[float, float], float],
        temperature_profile: Callable[[float, float], float],
        particle_count: int = 10_000,
        seed: int | None = None,
    ) -> None:
        rng = np.random.default_rng(seed)
        self.particles.clear()

        for _ in range(particle_count):
            # Sample cylindrical coordinates inside the torus volume
            r = self.minor_radius * rng.random() ** (1 / 3)
            theta = rng.uniform(0, 2 * np.pi)
            z = rng.uniform(-self.minor_radius, self.minor_radius)
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            position = np.array([x, y, z])

            density = density_profile(r, z)
            temperature = temperature_profile(r, z)
            if density <= 0 or temperature <= 0:
                continue

            ion_species = self.species["deuteron"] if rng.random() < 0.5 else self.species["triton"]
            electron_species = self.species["electron"]

            kT = temperature * 1.602e-19
            v_th_ion = math.sqrt(2 * kT / ion_species.mass)
            v_th_e = math.sqrt(2 * kT / electron_species.mass)
            ion_velocity = rng.normal(0, v_th_ion / math.sqrt(3), size=3)
            electron_velocity = rng.normal(0, v_th_e / math.sqrt(3), size=3)

            self.particles.append(Particle(ion_species, position.copy(), ion_velocity))
            self.particles.append(Particle(electron_species, position.copy(), electron_velocity))

    # ------------------------------------------------------------------
    # Core update loop
    # ------------------------------------------------------------------
    def run(self, controller: EnergyController | None = None) -> None:
        output_interval = 1e-4
        next_output = output_interval

        while self.time < self.t_max and not self.abort_triggered:
            self._step(controller)
            if self.time >= next_output:
                self.log.append(self._diagnostics())
                next_output += output_interval

    def _step(self, controller: EnergyController | None) -> None:
        self.field_solver.reset()
        for particle in self.particles:
            self.field_solver.deposit(particle)
        self.field_solver.solve()

        self._push_particles()
        self._handle_collisions()
        self._apply_losses()
        self._fusion_sampling()
        boundary = self._apply_boundaries()

        if controller is not None and self.particles:
            plasma_energy = sum(p.kinetic_energy() for p in self.particles)
            command = controller.compute_command(plasma_energy, self.time)
            self._apply_heating(command)

        self._adapt_timestep()
        self.time += self.dt

        if boundary.deposited_energy >= self.max_wall_energy:
            self.abort_triggered = True

    # ------------------------------------------------------------------
    def _push_particles(self) -> None:
        for particle in self.particles:
            electric, magnetic = self.field_solver.gather(particle)
            q_over_m = particle.species.charge / particle.species.mass

            v_minus = particle.velocity + q_over_m * electric * (0.5 * self.dt)
            t_vec = q_over_m * magnetic * (0.5 * self.dt)
            t_mag2 = np.dot(t_vec, t_vec)
            v_prime = v_minus + np.cross(v_minus, t_vec)
            s_vec = 2 * t_vec / (1 + t_mag2)
            v_plus = v_minus + np.cross(v_prime, s_vec)
            new_velocity = v_plus + q_over_m * electric * (0.5 * self.dt)
            new_velocity[2] -= G * self.dt
            particle.set_velocity(new_velocity)
            particle.advance(particle.velocity * self.dt)

    def _handle_collisions(self) -> None:
        rng = np.random.default_rng()
        angle_std = 0.05
        for particle in self.particles:
            speed = np.linalg.norm(particle.velocity)
            if speed == 0:
                continue
            axis = rng.normal(size=3)
            axis -= axis.dot(particle.velocity) * particle.velocity / (speed ** 2)
            axis_norm = np.linalg.norm(axis)
            if axis_norm == 0:
                continue
            axis /= axis_norm
            angle = rng.normal(0, angle_std)
            cos_a = math.cos(angle)
            sin_a = math.sin(angle)
            v_par = np.dot(particle.velocity, axis) * axis
            v_perp = particle.velocity - v_par
            rotated = v_par + cos_a * v_perp + sin_a * np.cross(axis, v_perp)
            particle.set_velocity(rotated)

    def _apply_losses(self) -> None:
        if self.energy_loss_rate <= 0:
            return
        damping = math.sqrt(max(0.0, 1.0 - self.energy_loss_rate * self.dt))
        for particle in self.particles:
            particle.set_velocity(particle.velocity * damping)

    def _apply_heating(self, power: float) -> None:
        if power <= 0 or not self.particles:
            return
        total_energy = sum(p.kinetic_energy() for p in self.particles)
        if total_energy <= 0:
            return
        fraction = power * self.dt / total_energy
        scale = math.sqrt(max(0.0, 1.0 + fraction))
        for particle in self.particles:
            particle.set_velocity(particle.velocity * scale)

    def _fusion_sampling(self) -> None:
        if self.fusion_reactivity <= 0:
            return
        rng = np.random.default_rng()
        new_particles: list[Particle] = []
        survivors: list[Particle] = []
        for particle in self.particles:
            if particle.species.name not in {"deuteron", "triton"}:
                survivors.append(particle)
                continue
            partner = rng.choice(self.particles)
            if partner.species.name == particle.species.name:
                survivors.append(particle)
                continue
            relative_velocity = np.linalg.norm(particle.velocity - partner.velocity)
            probability = self.fusion_reactivity * self.dt * max(relative_velocity, 0) / (self.grid_spacing[0] ** 3)
            if rng.random() < probability:
                alpha_species = self.species["alpha"]
                alpha_energy = 3.5e6 * 1.602e-19
                alpha_speed = math.sqrt(2 * alpha_energy / alpha_species.mass)
                direction = rng.normal(size=3)
                direction /= np.linalg.norm(direction)
                alpha_velocity = alpha_speed * direction
                new_particles.append(Particle(alpha_species, particle.position.copy(), alpha_velocity))
                self.total_neutron_energy += 14.1e6 * 1.602e-19
            else:
                survivors.append(particle)
        self.particles = survivors + new_particles

    def _apply_boundaries(self) -> BoundaryResult:
        result = BoundaryResult()
        survivors: list[Particle] = []
        for particle in self.particles:
            r = cylindrical_radius(particle.position, self.major_radius)
            z = particle.position[2]
            if r >= self.minor_radius or abs(z) >= self.minor_radius:
                energy = particle.kinetic_energy()
                result.deposited_energy += energy
                if z < -0.8 * self.minor_radius:
                    result.divertor_energy += energy
                else:
                    result.wall_energy += energy
                result.lost_particles += 1
            else:
                survivors.append(particle)
        self.cumulative_wall_energy += result.deposited_energy
        self.particles = survivors
        return result

    def _adapt_timestep(self) -> None:
        max_speed = max((np.linalg.norm(p.velocity) for p in self.particles), default=0.0)
        if max_speed > 0.5 * C and self.dt > 1e-9:
            self.dt /= 2
        elif max_speed < 0.1 * C and self.dt < 1e-7:
            self.dt *= 2

    # ------------------------------------------------------------------
    def _diagnostics(self) -> dict[str, float]:
        energy = sum(p.kinetic_energy() for p in self.particles)
        return {
            "time": self.time,
            "energy": energy,
            "particle_count": len(self.particles),
            "wall_energy": self.cumulative_wall_energy,
            "neutron_energy": self.total_neutron_energy,
        }


# ---------------------------------------------------------------------------
# Optional IMAS logging helper
# ---------------------------------------------------------------------------


class ImasLogger:
    """Persist simulation results to IMAS if the library is available."""

    def __init__(self, location: str) -> None:
        self.location = location
        self.enabled = _has_imas()
        if not self.enabled:
            self._fallback_records: list[dict[str, float]] = []
        else:
            import imas  # type: ignore
            self._db = imas.DBEntry(location, "w")
            self._db.open()

    def record(self, entry: dict[str, float]) -> None:
        if not self.enabled:
            self._fallback_records.append(entry)
            return
        from imas import ids  # type: ignore

        summary = ids(summary=True)
        summary.time = [entry["time"]]
        summary.global_quantities.fusion_energy = entry["neutron_energy"]
        summary.global_quantities.wall_energy = entry["wall_energy"]
        summary.global_quantities.plasma_duration = entry["time"]
        summary.global_quantities.max_plasma_energy = entry["energy"]
        self._db.put(summary)

    def close(self) -> None:
        if not self.enabled:
            return
        self._db.close()

    @property
    def fallback_records(self) -> list[dict[str, float]]:
        if self.enabled:
            raise RuntimeError("Fallback records only exist when IMAS is unavailable")
        return self._fallback_records
