"""Physical constants used throughout the tokamak simulation."""

import numpy as np

qe: float = 1.602176634e-19
"""Elementary charge (Coulombs)."""

me: float = 9.1093837015e-31
"""Electron mass (kg)."""

mp: float = 1.67262192369e-27
"""Proton mass (kg)."""

c: float = 299_792_458.0
"""Speed of light in vacuum (m/s)."""

MU0: float = 4e-7 * np.pi
"""Permeability of free space (H/m)."""

EPS0: float = 1.0 / (MU0 * c ** 2)
"""Permittivity of free space (F/m)."""

GRAVITY: float = 9.81
"""Gravitational acceleration (m/s^2)."""

DT_FUSION_ENERGY_J: float = 17.6e6 * qe
"""Energy released per D-T fusion reaction (Joules)."""

ALPHA_BIRTH_ENERGY_J: float = 3.5e6 * qe
"""Kinetic energy of the charged alpha product (Joules)."""

NEUTRON_ENERGY_J: float = DT_FUSION_ENERGY_J - ALPHA_BIRTH_ENERGY_J
"""Energy carried away by the neutron (Joules)."""
