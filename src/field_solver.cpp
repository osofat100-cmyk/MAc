#include "sim/field_solver.hpp"

#include <array>
#include <utility>

namespace sim {

FieldSolver::FieldSolver(const GridShape& shape, const GridSpacing& spacing, double toroidal_field)
    : shape_(shape), spacing_(spacing), toroidal_field_(toroidal_field) {}

void FieldSolver::reset() {
  // Placeholder implementation.
}

void FieldSolver::deposit(const Vec3& /*position*/, double /*charge*/, const Vec3& /*velocity*/) {
  // Placeholder implementation.
}

void FieldSolver::solve() {
  // Placeholder implementation.
}

std::pair<Vec3, Vec3> FieldSolver::gather(const Vec3& /*position*/) const {
  return {Vec3{0.0, 0.0, 0.0}, Vec3{0.0, 0.0, toroidal_field_}};
}

}  // namespace sim
