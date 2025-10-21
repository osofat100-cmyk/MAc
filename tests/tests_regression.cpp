#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include "sim/simulation.hpp"

namespace {

using sim::Vec3;

Vec3 JsonToVec3(const nlohmann::json& node) {
  return Vec3{
      node.at(0).get<double>(),
      node.at(1).get<double>(),
      node.at(2).get<double>(),
  };
}

void ExpectVecNear(const Vec3& actual, const Vec3& expected, double rel_tol = 1e-9, double abs_tol = 1e-12) {
  for (int i = 0; i < 3; ++i) {
    const double exp_val = expected[i];
    const double tolerance = std::max(abs_tol, std::abs(exp_val) * rel_tol);
    EXPECT_NEAR(exp_val, actual[i], tolerance);
  }
}

const sim::Species& SpeciesFromName(const std::string& name) {
  if (name == "electron") {
    return sim::electron;
  }
  if (name == "deuteron") {
    return sim::deuteron;
  }
  if (name == "triton") {
    return sim::triton;
  }
  if (name == "alpha") {
    return sim::alpha;
  }
  throw std::runtime_error("Unknown species name: " + name);
}

std::filesystem::path SnapshotPath(const std::string& filename) {
  static const std::filesystem::path base = std::filesystem::path{TEST_DATA_DIR};
  return base / filename;
}

}  // namespace

TEST(FieldSolverRegression, MatchesPythonSnapshot) {
  const std::filesystem::path path = SnapshotPath("python_field_snapshot.json");
  std::ifstream input(path);
  ASSERT_TRUE(input.is_open()) << "Unable to open snapshot at " << path.string();

  nlohmann::json snapshot;
  input >> snapshot;

  const auto& cfg = snapshot.at("config");
  sim::FieldSolver::GridShape shape{
      cfg.at("grid_shape").at(0).get<int>(),
      cfg.at("grid_shape").at(1).get<int>(),
      cfg.at("grid_shape").at(2).get<int>(),
  };
  sim::FieldSolver::GridSpacing spacing{
      cfg.at("grid_spacing").at(0).get<double>(),
      cfg.at("grid_spacing").at(1).get<double>(),
      cfg.at("grid_spacing").at(2).get<double>(),
  };
  const double toroidal_field = cfg.at("toroidal_field").get<double>();

  sim::FieldSolver solver(shape, spacing, toroidal_field);

  for (const auto& particle : snapshot.at("particles")) {
    const Vec3 position = JsonToVec3(particle.at("position"));
    const Vec3 velocity = JsonToVec3(particle.at("velocity"));
    const double charge = particle.at("charge").get<double>();
    solver.deposit(position, charge, velocity);
  }
  solver.solve();

  for (const auto& sample : snapshot.at("samples")) {
    const Vec3 position = JsonToVec3(sample.at("position"));
    const auto [electric, magnetic] = solver.gather(position);
    ExpectVecNear(electric, JsonToVec3(sample.at("electric")), 1e-6, 1e-12);
    ExpectVecNear(magnetic, JsonToVec3(sample.at("magnetic")), 1e-7, 1e-12);
  }
}

TEST(SimulationRegression, EnergyLoopsMatchPythonSnapshot) {
  const std::filesystem::path path = SnapshotPath("python_energy_snapshot.json");
  std::ifstream input(path);
  ASSERT_TRUE(input.is_open()) << "Unable to open snapshot at " << path.string();

  nlohmann::json snapshot;
  input >> snapshot;

  sim::SimulationConfig config;
  config.dt = snapshot.at("config").at("dt").get<double>();
  config.energy_loss_rate = snapshot.at("config").at("energy_loss_rate").get<double>();
  config.max_time = config.dt * 10.0;  // ensure we can heat without hitting completion.

  sim::Simulation simulation(config);

  const auto& particles = snapshot.at("particles");
  for (const auto& particle : particles) {
    const auto& species = SpeciesFromName(particle.at("species").get<std::string>());
    const Vec3 position = JsonToVec3(particle.at("position"));
    const Vec3 velocity = JsonToVec3(particle.at("initial_velocity"));
    const double weight = particle.at("weight").get<double>();
    simulation.add_particle(species, position, velocity, weight);
  }

  simulation.apply_losses();

  for (std::size_t i = 0; i < particles.size(); ++i) {
    const Vec3 expected = JsonToVec3(particles.at(i).at("after_losses"));
    ExpectVecNear(simulation.particles().velocity(i), expected, 1e-9, 1e-9);
  }

  const double power = snapshot.at("power").get<double>();
  simulation.apply_heating(power);

  for (std::size_t i = 0; i < particles.size(); ++i) {
    const Vec3 expected = JsonToVec3(particles.at(i).at("after_heating"));
    ExpectVecNear(simulation.particles().velocity(i), expected, 1e-9, 1e-6);
  }

  const double expected_energy = snapshot.at("total_energy").get<double>();
  EXPECT_NEAR(expected_energy, simulation.total_energy(), std::max(1e-12, expected_energy * 1e-6));
}
