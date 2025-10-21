#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "sim/constants.hpp"
#include "sim/species.hpp"

namespace sim {

using Vec3 = std::array<double, 3>;

class ParticleState {
 public:
  ParticleState() = default;

  [[nodiscard]] std::size_t size() const noexcept { return positions_.size(); }

  std::size_t emplace_particle(const Species& species, const Vec3& position, const Vec3& velocity, double weight = 1.0);

  [[nodiscard]] const Vec3& position(std::size_t index) const { return positions_[index]; }
  [[nodiscard]] const Vec3& velocity(std::size_t index) const { return velocities_[index]; }
  [[nodiscard]] double weight(std::size_t index) const { return weights_[index]; }
  [[nodiscard]] double gamma(std::size_t index) const { return gamma_[index]; }
  [[nodiscard]] const Species& species(std::size_t index) const { return *species_[index]; }

  void set_position(std::size_t index, const Vec3& value);
  void set_velocity(std::size_t index, const Vec3& value);
  void set_weight(std::size_t index, double weight) { weights_[index] = weight; }
  void set_species(std::size_t index, const Species& species) { species_[index] = &species; }

  void update_gamma(std::size_t index);
  void update_all_gamma();

  [[nodiscard]] double kinetic_energy(std::size_t index) const;
  [[nodiscard]] double total_kinetic_energy() const;

  void scale_velocity(std::size_t index, double factor);
  void scale_all_velocities(double factor);

 private:
  static double compute_gamma(const Vec3& velocity);

  std::vector<Vec3> positions_;
  std::vector<Vec3> velocities_;
  std::vector<double> weights_;
  std::vector<double> gamma_;
  std::vector<const Species*> species_;
};

}  // namespace sim
