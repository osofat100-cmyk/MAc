#pragma once

#include <array>
#include <optional>
#include <string>

#include "sim/field_solver.hpp"
#include "sim/particle.hpp"
#include "sim/species.hpp"

namespace sim {

struct SimulationConfig {
  double major_radius = 1.6;
  double minor_radius = 0.45;
  double toroidal_field = 3.5;
  double dt = 1e-8;
  double max_time = 5e-3;
  double energy_loss_rate = 0.5;
  double fusion_cross_section = 1e-21;
  double max_wall_energy = 5e7;
  std::array<int, 3> grid_shape{24, 24, 24};
  std::array<double, 3> grid_spacing{0.02, 0.02, 0.02};
};

class Simulation {
 public:
  explicit Simulation(const SimulationConfig& config = {});

  [[nodiscard]] const SimulationConfig& config() const noexcept { return config_; }
  [[nodiscard]] double time() const noexcept { return time_; }
  [[nodiscard]] double wall_energy() const noexcept { return wall_energy_; }
  [[nodiscard]] double fusion_energy() const noexcept { return fusion_energy_; }
  [[nodiscard]] const std::optional<std::string>& abort_reason() const noexcept { return abort_reason_; }

  ParticleState& particles() noexcept { return particles_; }
  const ParticleState& particles() const noexcept { return particles_; }

  void clear_particles();
  std::size_t add_particle(const Species& species, const Vec3& position, const Vec3& velocity, double weight = 1.0);

  void push_particles();
  double total_energy() const { return particles_.total_kinetic_energy(); }

 private:
  SimulationConfig config_;
  ParticleState particles_;
  FieldSolver field_solver_;
  double time_{0.0};
  double wall_energy_{0.0};
  double fusion_energy_{0.0};
  std::optional<std::string> abort_reason_;
};

}  // namespace sim
