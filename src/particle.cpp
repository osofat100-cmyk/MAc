#include "sim/particle.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace sim {

std::size_t ParticleState::emplace_particle(const Species& species, const Vec3& position, const Vec3& velocity, double weight) {
  positions_.push_back(position);
  velocities_.push_back(velocity);
  weights_.push_back(weight);
  species_.push_back(&species);
  gamma_.push_back(compute_gamma(velocity));
  return positions_.size() - 1;
}

void ParticleState::set_position(std::size_t index, const Vec3& value) {
  positions_[index] = value;
}

void ParticleState::set_velocity(std::size_t index, const Vec3& value) {
  velocities_[index] = value;
  gamma_[index] = compute_gamma(value);
}

void ParticleState::update_gamma(std::size_t index) {
  gamma_[index] = compute_gamma(velocities_[index]);
}

void ParticleState::update_all_gamma() {
  for (std::size_t i = 0; i < velocities_.size(); ++i) {
    gamma_[i] = compute_gamma(velocities_[i]);
  }
}

double ParticleState::kinetic_energy(std::size_t index) const {
  const auto* species = species_[index];
  return (gamma_[index] - 1.0) * species->mass * constants::c * constants::c * weights_[index];
}

double ParticleState::total_kinetic_energy() const {
  double total = 0.0;
  for (std::size_t i = 0; i < positions_.size(); ++i) {
    const auto* species = species_[i];
    total += (gamma_[i] - 1.0) * species->mass * constants::c * constants::c * weights_[i];
  }
  return total;
}

void ParticleState::scale_velocity(std::size_t index, double factor) {
  auto& velocity = velocities_[index];
  velocity[0] *= factor;
  velocity[1] *= factor;
  velocity[2] *= factor;
  gamma_[index] = compute_gamma(velocity);
}

void ParticleState::scale_all_velocities(double factor) {
  for (std::size_t i = 0; i < velocities_.size(); ++i) {
    auto& v = velocities_[i];
    v[0] *= factor;
    v[1] *= factor;
    v[2] *= factor;
    gamma_[i] = compute_gamma(v);
  }
}

double ParticleState::compute_gamma(const Vec3& velocity) {
  const double c2 = constants::c * constants::c;
  const double v2 = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
  const double clamped = std::min(v2, 0.999999 * c2);
  return 1.0 / std::sqrt(1.0 - clamped / c2);
}

}  // namespace sim
