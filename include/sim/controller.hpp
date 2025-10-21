#pragma once

#include <algorithm>

namespace sim {

class EnergyController {
 public:
  EnergyController(double target_energy, double critical_energy, double proportional_gain = 0.2,
                   double nonlinear_gain = 0.8, double min_power = 0.0, double max_power = 100e6)
      : target_energy_(target_energy),
        critical_energy_(critical_energy),
        proportional_gain_(proportional_gain),
        nonlinear_gain_(nonlinear_gain),
        min_power_(min_power),
        max_power_(max_power) {}

  [[nodiscard]] double command(double energy) const {
    const double error = target_energy_ - energy;
    const double nonlinear = energy * (1.0 - energy / critical_energy_);
    const double power = proportional_gain_ * error + nonlinear_gain_ * nonlinear;
    return std::clamp(power, min_power_, max_power_);
  }

  void set_gains(double proportional, double nonlinear) {
    proportional_gain_ = proportional;
    nonlinear_gain_ = nonlinear;
  }

  void set_limits(double min_power, double max_power) {
    min_power_ = min_power;
    max_power_ = max_power;
  }

 private:
  double target_energy_;
  double critical_energy_;
  double proportional_gain_;
  double nonlinear_gain_;
  double min_power_;
  double max_power_;
};

}  // namespace sim
