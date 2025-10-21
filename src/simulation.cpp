#include "sim/simulation.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
#include <vector>

#include "sim/constants.hpp"

namespace sim {

namespace {

[[nodiscard]] Vec3 add(const Vec3& a, const Vec3& b) {
  return Vec3{a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

[[nodiscard]] Vec3 subtract(const Vec3& a, const Vec3& b) {
  return Vec3{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
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

[[nodiscard]] double norm(const Vec3& v) {
  return std::sqrt(dot(v, v));
}

[[nodiscard]] Vec3 normalise(const Vec3& v) {
  const double n = norm(v);
  if (n == 0.0) {
    return Vec3{0.0, 0.0, 0.0};
  }
  return scale(v, 1.0 / n);
}

[[nodiscard]] double clamp_density(double density) {
  return density <= 0.0 ? 0.0 : density;
}

}  // namespace

Simulation::Simulation(const SimulationConfig& config)
    : config_(config),
      field_solver_(config.grid_shape, config.grid_spacing, config.toroidal_field),
      rng_(config.rng_seed != 0 ? config.rng_seed : std::random_device{}()) {}

void Simulation::clear_particles() {
  particles_ = ParticleState{};
}

std::size_t Simulation::add_particle(const Species& species, const Vec3& position, const Vec3& velocity,
                                     double weight) {
  return particles_.emplace_particle(species, position, velocity, weight);
}

void Simulation::initialise_plasma(const Profile& density_profile, const Profile& temperature_profile,
                                   std::size_t macro_particles) {
  initialise_plasma(density_profile, temperature_profile, macro_particles, rng_);
}

void Simulation::initialise_plasma(const Profile& density_profile, const Profile& temperature_profile,
                                   std::size_t macro_particles, std::mt19937_64& rng) {
  clear_particles();

  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  std::normal_distribution<double> normal01(0.0, 1.0);

  const double volume = 2.0 * std::numbers::pi * config_.major_radius * std::numbers::pi *
                        config_.minor_radius * config_.minor_radius;
  const double peak_density = std::max(density_profile(0.0, 0.0), 1e16);

  while (particles_.size() < macro_particles) {
    const double u_r = uniform01(rng);
    const double r = config_.minor_radius * std::pow(u_r, 1.0 / 3.0);
    const double theta = 2.0 * std::numbers::pi * uniform01(rng);
    const double z = (uniform01(rng) * 2.0 - 1.0) * config_.minor_radius;

    const double x = r * std::cos(theta);
    const double y = r * std::sin(theta);

    const double density = clamp_density(density_profile(r, z));
    if (density <= 0.0) {
      continue;
    }
    if (uniform01(rng) > density / peak_density) {
      continue;
    }

    const double temperature = std::max(temperature_profile(r, z), 1.0);
    const Species& ion_species = (uniform01(rng) < 0.5) ? deuteron : triton;
    const double ion_vth = std::sqrt(2.0 * temperature * constants::qe / ion_species.mass);
    const double electron_vth = std::sqrt(2.0 * temperature * constants::qe / electron.mass);

    const double sigma_ion = ion_vth / std::sqrt(3.0);
    const double sigma_electron = electron_vth / std::sqrt(3.0);

    Vec3 ion_velocity{sigma_ion * normal01(rng), sigma_ion * normal01(rng), sigma_ion * normal01(rng)};
    Vec3 electron_velocity{
        sigma_electron * normal01(rng),
        sigma_electron * normal01(rng),
        sigma_electron * normal01(rng),
    };

    Vec3 position{x, y, z};

    particles_.emplace_particle(ion_species, position, ion_velocity);
    particles_.emplace_particle(electron, position, electron_velocity);
  }

  const std::size_t count = particles_.size();
  if (count == 0) {
    return;
  }

  const double weight = density_profile(0.0, 0.0) * volume / static_cast<double>(count);
  for (std::size_t i = 0; i < count; ++i) {
    particles_.set_weight(i, weight);
  }
}

void Simulation::deposit_fields() {
  field_solver_.reset();
  const std::size_t count = particles_.size();
  for (std::size_t i = 0; i < count; ++i) {
    const Vec3& position = particles_.position(i);
    const Vec3& velocity = particles_.velocity(i);
    const Species& species = particles_.species(i);
    const double charge = species.charge * particles_.weight(i);
    field_solver_.deposit(position, charge, velocity);
  }
  field_solver_.solve();
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

void Simulation::handle_collisions() {
  std::normal_distribution<double> normal_axis(0.0, 1.0);
  std::normal_distribution<double> normal_angle(0.0, 0.05);

  const std::size_t count = particles_.size();
  for (std::size_t i = 0; i < count; ++i) {
    Vec3 velocity = particles_.velocity(i);
    const double speed2 = dot(velocity, velocity);
    if (speed2 == 0.0) {
      continue;
    }

    Vec3 axis{normal_axis(rng_), normal_axis(rng_), normal_axis(rng_)};
    const double proj = dot(axis, velocity);
    axis = subtract(axis, scale(velocity, proj / speed2));

    const double axis_norm = norm(axis);
    if (axis_norm == 0.0) {
      continue;
    }
    axis = scale(axis, 1.0 / axis_norm);

    const double angle = normal_angle(rng_);
    const double cos_a = std::cos(angle);
    const double sin_a = std::sin(angle);

    Vec3 rotated = add(
        add(scale(velocity, cos_a), scale(cross(axis, velocity), sin_a)),
        scale(axis, dot(axis, velocity) * (1.0 - cos_a)));

    particles_.set_velocity(i, rotated);
  }
}

void Simulation::apply_losses() {
  if (config_.energy_loss_rate <= 0.0) {
    return;
  }
  const double factor = std::exp(-config_.energy_loss_rate * config_.dt / 2.0);
  if (!std::isfinite(factor) || factor <= 0.0) {
    return;
  }
  particles_.scale_all_velocities(factor);
}

void Simulation::fusion_events() {
  const std::size_t count = particles_.size();
  if (count < 2) {
    return;
  }

  const double cell_volume = static_cast<double>(config_.grid_shape[0]) * config_.grid_spacing[0] *
                             static_cast<double>(config_.grid_shape[1]) * config_.grid_spacing[1] *
                             static_cast<double>(config_.grid_shape[2]) * config_.grid_spacing[2];
  const double rate = (config_.fusion_cross_section * config_.dt) / std::max(cell_volume, 1e-12);

  std::uniform_real_distribution<double> uniform01(0.0, 1.0);
  std::normal_distribution<double> normal01(0.0, 1.0);

  std::vector<char> removed(count, 0);
  struct AlphaCandidate {
    Vec3 position;
    Vec3 velocity;
    double weight;
  };
  std::vector<AlphaCandidate> alphas;
  double neutron_energy = 0.0;

  for (std::size_t i = 0; i < count; ++i) {
    if (removed[i]) {
      continue;
    }
    const Species& si = particles_.species(i);
    if (&si != &deuteron && &si != &triton) {
      continue;
    }
    for (std::size_t j = i + 1; j < count; ++j) {
      if (removed[j]) {
        continue;
      }
      const Species& sj = particles_.species(j);
      const bool is_pair = ((&si == &deuteron && &sj == &triton) || (&si == &triton && &sj == &deuteron));
      if (!is_pair) {
        continue;
      }

      if (uniform01(rng_) >= rate) {
        continue;
      }

      Vec3 direction{normal01(rng_), normal01(rng_), normal01(rng_)};
      direction = normalise(direction);
      if (norm(direction) == 0.0) {
        continue;
      }

      const double speed = std::sqrt(2.0 * constants::alpha_birth_energy_j / alpha.mass);
      Vec3 velocity = scale(direction, speed);

      Vec3 position = scale(add(particles_.position(i), particles_.position(j)), 0.5);
      const double weight = 0.5 * (particles_.weight(i) + particles_.weight(j));

      alphas.push_back(AlphaCandidate{position, velocity, weight});
      neutron_energy += constants::neutron_energy_j;
      removed[i] = 1;
      removed[j] = 1;
      break;
    }
  }

  if (alphas.empty() && std::none_of(removed.begin(), removed.end(), [](char flag) { return flag != 0; })) {
    return;
  }

  ParticleState survivors;
  survivors.reserve(count - static_cast<std::size_t>(std::count(removed.begin(), removed.end(), 1)) + alphas.size());

  for (std::size_t i = 0; i < count; ++i) {
    if (removed[i]) {
      continue;
    }
    const Species& species = particles_.species(i);
    survivors.emplace_particle(species, particles_.position(i), particles_.velocity(i), particles_.weight(i));
  }

  for (const auto& alpha_candidate : alphas) {
    survivors.emplace_particle(alpha, alpha_candidate.position, alpha_candidate.velocity, alpha_candidate.weight);
  }

  particles_ = std::move(survivors);
  fusion_energy_ += neutron_energy;
}

void Simulation::apply_boundary() {
  const std::size_t count = particles_.size();
  if (count == 0) {
    return;
  }

  ParticleState survivors;
  survivors.reserve(count);

  double wall_energy_gain = 0.0;
  const double minor_radius = config_.minor_radius;
  for (std::size_t i = 0; i < count; ++i) {
    const Vec3& position = particles_.position(i);
    const double radial = std::sqrt(position[0] * position[0] + position[1] * position[1]);
    if (radial >= minor_radius || std::abs(position[2]) >= minor_radius) {
      wall_energy_gain += particles_.kinetic_energy(i);
      continue;
    }
    const Species& species = particles_.species(i);
    survivors.emplace_particle(species, position, particles_.velocity(i), particles_.weight(i));
  }

  particles_ = std::move(survivors);
  wall_energy_ += wall_energy_gain;

  if (wall_energy_ > config_.max_wall_energy) {
    abort_reason_ = "wall energy limit exceeded";
  }
}

void Simulation::apply_heating(double power) {
  if (power <= 0.0 || particles_.size() == 0) {
    return;
  }
  const double total = total_energy();
  if (total <= 0.0) {
    return;
  }
  const double scale_factor = std::sqrt(1.0 + power * config_.dt / total);
  if (!std::isfinite(scale_factor) || scale_factor <= 0.0) {
    return;
  }
  particles_.scale_all_velocities(scale_factor);
}

void Simulation::step() {
  if (abort_reason_) {
    return;
  }

  deposit_fields();
  push_particles();
  handle_collisions();
  apply_losses();
  fusion_events();
  apply_boundary();

  time_ += config_.dt;
  if (time_ >= config_.max_time) {
    abort_reason_ = "completed";
  }
}

void Simulation::run(const EnergyController* controller) {
  energy_history_.clear();
  time_history_.clear();

  while (!abort_reason_) {
    step();
    const double energy = total_energy();
    energy_history_.push_back(energy);
    time_history_.push_back(time_);
    if (controller) {
      const double command = controller->command(energy);
      apply_heating(command);
    }
  }
}

}  // namespace sim
