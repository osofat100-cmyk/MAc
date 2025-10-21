#pragma once

namespace sim::constants {

inline constexpr double qe = 1.602176634e-19;
inline constexpr double me = 9.1093837015e-31;
inline constexpr double mp = 1.67262192369e-27;
inline constexpr double c = 299'792'458.0;
inline constexpr double mu0 = 4e-7 * 3.14159265358979323846;
inline constexpr double eps0 = 1.0 / (mu0 * c * c);
inline constexpr double gravity = 9.81;
inline constexpr double dt_fusion_energy_j = 17.6e6 * qe;
inline constexpr double alpha_birth_energy_j = 3.5e6 * qe;
inline constexpr double neutron_energy_j = dt_fusion_energy_j - alpha_birth_energy_j;

}  // namespace sim::constants
