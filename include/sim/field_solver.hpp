#pragma once

#include <array>

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
  GridShape shape_{};
  GridSpacing spacing_{};
  double toroidal_field_{};
};

}  // namespace sim
