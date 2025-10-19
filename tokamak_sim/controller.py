"""Control algorithms for regulating plasma energy."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class EnergyController:
    """Nonlinear energy confinement stabilizer.

    The controller implements the feedback law

    u(t) = Kp * (E_star - E) + alpha * E * (1 - E / E_c)

    Parameters
    ----------
    E_star:
        Desired stored energy setpoint in joules.
    E_crit:
        Critical stored energy limit in joules.
    Kp:
        Linear proportional gain.
    alpha:
        Nonlinear gain shaping the saturation characteristic.
    Gamma:
        Open-loop loss rate used for documentation and post-processing.
    max_power:
        Upper bound on commanded auxiliary heating (W).
    min_power:
        Lower bound on commanded heating. Negative power is clamped to zero by
        default because the model does not include active power extraction.
    """

    E_star: float
    E_crit: float
    Kp: float = 0.2
    alpha: float = 0.8
    Gamma: float = 0.5
    max_power: float = 1e8
    min_power: float = 0.0

    def compute_command(self, energy: float, time: float) -> float:
        """Return the requested heating power in watts."""

        command = self.Kp * (self.E_star - energy) + self.alpha * energy * (1.0 - energy / self.E_crit)
        if command < self.min_power:
            return self.min_power
        if command > self.max_power:
            return self.max_power
        return command
