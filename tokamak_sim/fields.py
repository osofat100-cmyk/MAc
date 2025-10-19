"""Simplified electromagnetic field solver for the tokamak simulation."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .constants import EPS0, MU0


@dataclass
class FieldSolver:
    """Very small 3-D quasi-static field solver."""

    grid_shape: tuple[int, int, int]
    grid_spacing: tuple[float, float, float]
    background_toroidal_field: float

    def __post_init__(self) -> None:
        self.Nx, self.Ny, self.Nz = self.grid_shape
        self.dx, self.dy, self.dz = self.grid_spacing
        self.rho = np.zeros(self.grid_shape, dtype=float)
        self.J = np.zeros((*self.grid_shape, 3), dtype=float)
        self.E = np.zeros((*self.grid_shape, 3), dtype=float)
        self.B = np.zeros((*self.grid_shape, 3), dtype=float)

    def reset(self) -> None:
        self.rho.fill(0.0)
        self.J.fill(0.0)

    def deposit(self, position: np.ndarray, charge: float, velocity: np.ndarray) -> None:
        ix = int(np.floor(position[0] / self.dx + self.Nx / 2))
        iy = int(np.floor(position[1] / self.dy + self.Ny / 2))
        iz = int(np.floor(position[2] / self.dz + self.Nz / 2))
        if 0 <= ix < self.Nx and 0 <= iy < self.Ny and 0 <= iz < self.Nz:
            self.rho[ix, iy, iz] += charge
            self.J[ix, iy, iz] += charge * velocity

    def solve(self) -> None:
        # Poisson solve using Jacobi iterations
        phi = np.zeros_like(self.rho)
        laplace_coeff = 1.0 / (2.0 / self.dx ** 2 + 2.0 / self.dy ** 2 + 2.0 / self.dz ** 2)
        for _ in range(20):
            phi[1:-1, 1:-1, 1:-1] = (
                (phi[2:, 1:-1, 1:-1] + phi[:-2, 1:-1, 1:-1]) / self.dx ** 2
                + (phi[1:-1, 2:, 1:-1] + phi[1:-1, :-2, 1:-1]) / self.dy ** 2
                + (phi[1:-1, 1:-1, 2:] + phi[1:-1, 1:-1, :-2]) / self.dz ** 2
                + self.rho[1:-1, 1:-1, 1:-1] / EPS0
            ) * laplace_coeff
        Ex = (np.roll(phi, -1, axis=0) - np.roll(phi, 1, axis=0)) / (2 * self.dx)
        Ey = (np.roll(phi, -1, axis=1) - np.roll(phi, 1, axis=1)) / (2 * self.dy)
        Ez = (np.roll(phi, -1, axis=2) - np.roll(phi, 1, axis=2)) / (2 * self.dz)
        self.E[..., 0] = -Ex
        self.E[..., 1] = -Ey
        self.E[..., 2] = -Ez

        curl_x = (
            (np.roll(self.J[..., 2], -1, axis=1) - np.roll(self.J[..., 2], 1, axis=1)) / (2 * self.dy)
            - (np.roll(self.J[..., 1], -1, axis=2) - np.roll(self.J[..., 1], 1, axis=2)) / (2 * self.dz)
        )
        curl_y = (
            (np.roll(self.J[..., 0], -1, axis=2) - np.roll(self.J[..., 0], 1, axis=2)) / (2 * self.dz)
            - (np.roll(self.J[..., 2], -1, axis=0) - np.roll(self.J[..., 2], 1, axis=0)) / (2 * self.dx)
        )
        curl_z = (
            (np.roll(self.J[..., 1], -1, axis=0) - np.roll(self.J[..., 1], 1, axis=0)) / (2 * self.dx)
            - (np.roll(self.J[..., 0], -1, axis=1) - np.roll(self.J[..., 0], 1, axis=1)) / (2 * self.dy)
        )
        self.B[..., 0] = MU0 * curl_x * self.dx
        self.B[..., 1] = MU0 * curl_y * self.dy + self.background_toroidal_field
        self.B[..., 2] = MU0 * curl_z * self.dz

    def gather(self, position: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ix = int(np.floor(position[0] / self.dx + self.Nx / 2))
        iy = int(np.floor(position[1] / self.dy + self.Ny / 2))
        iz = int(np.floor(position[2] / self.dz + self.Nz / 2))
        if 0 <= ix < self.Nx and 0 <= iy < self.Ny and 0 <= iz < self.Nz:
            return self.E[ix, iy, iz].copy(), self.B[ix, iy, iz].copy()
        return np.zeros(3), np.zeros(3)
