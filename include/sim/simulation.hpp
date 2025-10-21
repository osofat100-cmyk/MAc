#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "sim/controller.hpp"
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
  std::uint64_t rng_seed = 0;
};

class Simulation {
 public:
  explicit Simulation(const SimulationConfig& config = {});

  using Profile = std::function<double(double, double)>;

  [[nodiscard]] const SimulationConfig& config() const noexcept { return config_; }
  [[nodiscard]] double time() const noexcept { return time_; }
  [[nodiscard]] double wall_energy() const noexcept { return wall_energy_; }
  [[nodiscard]] double fusion_energy() const noexcept { return fusion_energy_; }
  [[nodiscard]] const std::optional<std::string>& abort_reason() const noexcept { return abort_reason_; }

  ParticleState& particles() noexcept { return particles_; }
  const ParticleState& particles() const noexcept { return particles_; }

  void clear_particles();
  std::size_t add_particle(const Species& species, const Vec3& position, const Vec3& velocity,
                           double weight = 1.0);

  void initialise_plasma(const Profile& density_profile, const Profile& temperature_profile,
                         std::size_t macro_particles = 5000);
  void initialise_plasma(const Profile& density_profile, const Profile& temperature_profile,
                         std::size_t macro_particles, std::mt19937_64& rng);

  void deposit_fields();
  void push_particles();
  void handle_collisions();
  void apply_losses();
  void fusion_events();
  void apply_boundary();
  void apply_heating(double power);
  void step();
  void run(const EnergyController* controller = nullptr);

  double total_energy() const { return particles_.total_kinetic_energy(); }

  [[nodiscard]] const std::vector<double>& energy_history() const noexcept { return energy_history_; }
  [[nodiscard]] const std::vector<double>& time_history() const noexcept { return time_history_; }

  [[nodiscard]] std::mt19937_64& rng() noexcept { return rng_; }
  [[nodiscard]] const std::mt19937_64& rng() const noexcept { return rng_; }

 private:
  SimulationConfig config_;
  ParticleState particles_;
  FieldSolver field_solver_;
  double time_{0.0};
  double wall_energy_{0.0};
  double fusion_energy_{0.0};
  std::optional<std::string> abort_reason_;
  std::mt19937_64 rng_;
  std::vector<double> energy_history_;
  std::vector<double> time_history_;
};

}  // namespace sim
