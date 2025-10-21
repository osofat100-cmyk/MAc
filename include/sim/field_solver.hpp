#pragma once

#include <array>
#include <utility>
#include <vector>

#include "sim/particle.hpp"

namespace sim {

class FieldSolver {
 public:
  using GridShape = std::array<int, 3>;
  using GridSpacing = std::array<double, 3>;

  FieldSolver(const GridShape& shape, const GridSpacing& spacing, double toroidal_field);

  void reset();
  void deposit(const Vec3& position, double charge, const Vec3& velocity);
  void solve();

  [[nodiscard]] std::pair<Vec3, Vec3> gather(const Vec3& position) const;

 private:
  [[nodiscard]] std::size_t index(int i, int j, int k) const noexcept;
  [[nodiscard]] static int wrap_index(int value, int extent) noexcept;
  [[nodiscard]] bool in_bounds(int i, int j, int k) const noexcept;
  void compute_electric_field(const std::vector<double>& potential);
  void compute_magnetic_field();

  GridShape shape_{};
  GridSpacing spacing_{};
  double toroidal_field_{};
  std::size_t total_size_{0};

  double dx_{1.0};
  double dy_{1.0};
  double dz_{1.0};
  double inv_two_dx_{0.0};
  double inv_two_dy_{0.0};
  double inv_two_dz_{0.0};

  std::vector<double> charge_density_;
  std::vector<Vec3> current_density_;
  std::vector<Vec3> electric_field_;
  std::vector<Vec3> magnetic_field_;
};

}  // namespace sim
