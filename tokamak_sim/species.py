"""Species definitions for the tokamak simulation."""

from dataclasses import dataclass

from .constants import qe, me, mp


@dataclass(frozen=True)
class Species:
    """Physical definition of a particle species."""

    name: str
    mass: float
    charge: float

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.name


ELECTRON = Species("electron", me, -qe)
DEUTERON = Species("deuteron", 2.0 * mp, +qe)
TRITON = Species("triton", 3.0 * mp, +qe)
ALPHA = Species("alpha", 4.0 * mp, 2.0 * qe)
