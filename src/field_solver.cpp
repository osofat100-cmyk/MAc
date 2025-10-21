#include "sim/field_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "sim/constants.hpp"

namespace sim {

namespace {

constexpr int kJacobiIterations = 20;

}  // namespace

FieldSolver::FieldSolver(const GridShape& shape, const GridSpacing& spacing, double toroidal_field)
    : shape_(shape), spacing_(spacing), toroidal_field_(toroidal_field) {
  total_size_ = static_cast<std::size_t>(shape_[0]) * static_cast<std::size_t>(shape_[1]) *
                static_cast<std::size_t>(shape_[2]);
  dx_ = spacing_[0];
  dy_ = spacing_[1];
  dz_ = spacing_[2];
  inv_two_dx_ = (dx_ != 0.0) ? 1.0 / (2.0 * dx_) : 0.0;
  inv_two_dy_ = (dy_ != 0.0) ? 1.0 / (2.0 * dy_) : 0.0;
  inv_two_dz_ = (dz_ != 0.0) ? 1.0 / (2.0 * dz_) : 0.0;

  charge_density_.assign(total_size_, 0.0);
  current_density_.assign(total_size_, Vec3{0.0, 0.0, 0.0});
  electric_field_.assign(total_size_, Vec3{0.0, 0.0, 0.0});
  magnetic_field_.assign(total_size_, Vec3{0.0, toroidal_field_, 0.0});
}

void FieldSolver::reset() {
  std::fill(charge_density_.begin(), charge_density_.end(), 0.0);
  std::fill(current_density_.begin(), current_density_.end(), Vec3{0.0, 0.0, 0.0});
  std::fill(electric_field_.begin(), electric_field_.end(), Vec3{0.0, 0.0, 0.0});
  std::fill(magnetic_field_.begin(), magnetic_field_.end(), Vec3{0.0, toroidal_field_, 0.0});
}

void FieldSolver::deposit(const Vec3& position, double charge, const Vec3& velocity) {
  int ix = static_cast<int>(std::floor(position[0] / dx_ + shape_[0] * 0.5));
  int iy = static_cast<int>(std::floor(position[1] / dy_ + shape_[1] * 0.5));
  int iz = static_cast<int>(std::floor(position[2] / dz_ + shape_[2] * 0.5));

  if (!in_bounds(ix, iy, iz)) {
    return;
  }

  const std::size_t idx = index(ix, iy, iz);
  charge_density_[idx] += charge;
  current_density_[idx][0] += charge * velocity[0];
  current_density_[idx][1] += charge * velocity[1];
  current_density_[idx][2] += charge * velocity[2];
}

void FieldSolver::solve() {
  if (total_size_ == 0) {
    return;
  }

  std::vector<double> phi(total_size_, 0.0);
  std::vector<double> phi_next(total_size_, 0.0);

  const double inv_dx2 = (dx_ != 0.0) ? 1.0 / (dx_ * dx_) : 0.0;
  const double inv_dy2 = (dy_ != 0.0) ? 1.0 / (dy_ * dy_) : 0.0;
  const double inv_dz2 = (dz_ != 0.0) ? 1.0 / (dz_ * dz_) : 0.0;
  const double laplace_coeff = 1.0 / (2.0 * inv_dx2 + 2.0 * inv_dy2 + 2.0 * inv_dz2 + 1e-16);
  const double rho_scale = 1.0 / constants::eps0;

  const int nx = shape_[0];
  const int ny = shape_[1];
  const int nz = shape_[2];

  for (int iter = 0; iter < kJacobiIterations; ++iter) {
    for (int i = 1; i < nx - 1; ++i) {
      for (int j = 1; j < ny - 1; ++j) {
        for (int k = 1; k < nz - 1; ++k) {
          const std::size_t idx = index(i, j, k);
          const double neighbors =
              (phi[index(i + 1, j, k)] + phi[index(i - 1, j, k)]) * inv_dx2 +
              (phi[index(i, j + 1, k)] + phi[index(i, j - 1, k)]) * inv_dy2 +
              (phi[index(i, j, k + 1)] + phi[index(i, j, k - 1)]) * inv_dz2 +
              charge_density_[idx] * rho_scale;
          phi_next[idx] = neighbors * laplace_coeff;
        }
      }
    }
    std::swap(phi, phi_next);
  }

  compute_electric_field(phi);
  compute_magnetic_field();
}

std::pair<Vec3, Vec3> FieldSolver::gather(const Vec3& position) const {
  int ix = static_cast<int>(std::floor(position[0] / dx_ + shape_[0] * 0.5));
  int iy = static_cast<int>(std::floor(position[1] / dy_ + shape_[1] * 0.5));
  int iz = static_cast<int>(std::floor(position[2] / dz_ + shape_[2] * 0.5));

  if (!in_bounds(ix, iy, iz)) {
    return {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, 0.0}};
  }

  const std::size_t idx = index(ix, iy, iz);
  return {electric_field_[idx], magnetic_field_[idx]};
}

std::size_t FieldSolver::index(int i, int j, int k) const noexcept {
  const int ny = shape_[1];
  const int nz = shape_[2];
  return static_cast<std::size_t>((i * ny + j) * nz + k);
}

int FieldSolver::wrap_index(int value, int extent) noexcept {
  const int mod = value % extent;
  return mod < 0 ? mod + extent : mod;
}

bool FieldSolver::in_bounds(int i, int j, int k) const noexcept {
  return (i >= 0 && i < shape_[0] && j >= 0 && j < shape_[1] && k >= 0 && k < shape_[2]);
}

void FieldSolver::compute_electric_field(const std::vector<double>& potential) {
  const int nx = shape_[0];
  const int ny = shape_[1];
  const int nz = shape_[2];

  for (int i = 0; i < nx; ++i) {
    const int ip = wrap_index(i + 1, nx);
    const int im = wrap_index(i - 1, nx);
    for (int j = 0; j < ny; ++j) {
      const int jp = wrap_index(j + 1, ny);
      const int jm = wrap_index(j - 1, ny);
      for (int k = 0; k < nz; ++k) {
        const int kp = wrap_index(k + 1, nz);
        const int km = wrap_index(k - 1, nz);
        const std::size_t idx = index(i, j, k);
        const double ex = -(potential[index(ip, j, k)] - potential[index(im, j, k)]) * inv_two_dx_;
        const double ey = -(potential[index(i, jp, k)] - potential[index(i, jm, k)]) * inv_two_dy_;
        const double ez = -(potential[index(i, j, kp)] - potential[index(i, j, km)]) * inv_two_dz_;
        electric_field_[idx] = Vec3{ex, ey, ez};
      }
    }
  }
}

void FieldSolver::compute_magnetic_field() {
  const int nx = shape_[0];
  const int ny = shape_[1];
  const int nz = shape_[2];

  for (int i = 0; i < nx; ++i) {
    const int ip = wrap_index(i + 1, nx);
    const int im = wrap_index(i - 1, nx);
    for (int j = 0; j < ny; ++j) {
      const int jp = wrap_index(j + 1, ny);
      const int jm = wrap_index(j - 1, ny);
      for (int k = 0; k < nz; ++k) {
        const int kp = wrap_index(k + 1, nz);
        const int km = wrap_index(k - 1, nz);
        const std::size_t idx = index(i, j, k);

        const Vec3& J_jp = current_density_[index(i, jp, k)];
        const Vec3& J_jm = current_density_[index(i, jm, k)];
        const Vec3& J_kp = current_density_[index(i, j, kp)];
        const Vec3& J_km = current_density_[index(i, j, km)];
        const Vec3& J_ip = current_density_[index(ip, j, k)];
        const Vec3& J_im = current_density_[index(im, j, k)];

        const double curl_x = (J_jp[2] - J_jm[2]) * inv_two_dy_ - (J_kp[1] - J_km[1]) * inv_two_dz_;
        const double curl_y = (J_kp[0] - J_km[0]) * inv_two_dz_ - (J_ip[2] - J_im[2]) * inv_two_dx_;
        const double curl_z = (J_ip[1] - J_im[1]) * inv_two_dx_ - (J_jp[0] - J_jm[0]) * inv_two_dy_;

        Vec3 B{
            constants::mu0 * curl_x * dx_,
            constants::mu0 * curl_y * dy_ + toroidal_field_,
            constants::mu0 * curl_z * dz_,
        };

        magnetic_field_[idx] = B;
      }
    }
  }
}

}  // namespace sim
