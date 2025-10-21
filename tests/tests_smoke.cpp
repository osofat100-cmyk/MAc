#include <cmath>
#include <gtest/gtest.h>

#include "sim/simulation.hpp"
#include "sim/species.hpp"

TEST(SimulationSmokeTest, BorisPushKeepsStateFinite) {
  sim::Simulation simulation;
  simulation.add_particle(sim::deuteron, sim::Vec3{0.0, 0.0, 0.0}, sim::Vec3{1e3, 0.0, 0.0});
  const double initial_energy = simulation.total_energy();
  ASSERT_GT(initial_energy, 0.0);

  EXPECT_NO_THROW(simulation.push_particles());

  const double updated_energy = simulation.total_energy();
  EXPECT_TRUE(std::isfinite(updated_energy));
  EXPECT_GT(updated_energy, 0.0);
}
