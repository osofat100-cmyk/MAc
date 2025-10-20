"""Non-linear plasma energy controller."""

from dataclasses import dataclass


@dataclass
class EnergyController:
    """Implements the nonlinear energy stabiliser law."""

    target_energy: float
    critical_energy: float
    proportional_gain: float = 0.2
    nonlinear_gain: float = 0.8
    min_power: float = 0.0
    max_power: float = 100e6

    def command(self, energy: float) -> float:
        error = self.target_energy - energy
        nonlinear = energy * (1.0 - energy / self.critical_energy)
        power = self.proportional_gain * error + self.nonlinear_gain * nonlinear
        return min(self.max_power, max(self.min_power, power))
