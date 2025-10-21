#include "sim/simulation.hpp"

#include <cmath>

#include "sim/constants.hpp"

namespace sim {

namespace {

[[nodiscard]] Vec3 add(const Vec3& a, const Vec3& b) {
  return Vec3{a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

[[nodiscard]] Vec3 scale(const Vec3& v, double factor) {
  return Vec3{v[0] * factor, v[1] * factor, v[2] * factor};
}

[[nodiscard]] Vec3 cross(const Vec3& a, const Vec3& b) {
  return Vec3{
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
  };
}

[[nodiscard]] double dot(const Vec3& a, const Vec3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

}  // namespace

Simulation::Simulation(const SimulationConfig& config)
    : config_(config), field_solver_(config.grid_shape, config.grid_spacing, config.toroidal_field) {}

void Simulation::clear_particles() {
  particles_ = ParticleState{};
}

std::size_t Simulation::add_particle(const Species& species, const Vec3& position, const Vec3& velocity, double weight) {
  return particles_.emplace_particle(species, position, velocity, weight);
}

void Simulation::push_particles() {
  const double dt = config_.dt;
  const double half_dt = 0.5 * dt;
  const std::size_t count = particles_.size();

  for (std::size_t i = 0; i < count; ++i) {
    const Species& species = particles_.species(i);
    Vec3 position = particles_.position(i);
    Vec3 velocity = particles_.velocity(i);

    const auto [E, B] = field_solver_.gather(position);
    const double q_over_m = species.charge / species.mass;

    Vec3 v_minus = add(velocity, scale(E, half_dt * q_over_m));
    Vec3 t = scale(B, half_dt * q_over_m);
    const double t_mag2 = dot(t, t);
    Vec3 v_prime = add(v_minus, cross(v_minus, t));
    Vec3 s = scale(t, 2.0 / (1.0 + t_mag2));
    Vec3 v_plus = add(v_minus, cross(v_prime, s));
    Vec3 new_velocity = add(v_plus, scale(E, half_dt * q_over_m));
    new_velocity[2] -= constants::gravity * dt;
    Vec3 new_position = add(position, scale(new_velocity, dt));

    particles_.set_velocity(i, new_velocity);
    particles_.set_position(i, new_position);
  }
}

}  // namespace sim
