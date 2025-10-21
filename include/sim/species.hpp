#pragma once

#include <string_view>

#include "sim/constants.hpp"

namespace sim {

struct Species {
  std::string_view name;
  double mass;
  double charge;
};

inline constexpr Species electron{"electron", constants::me, -constants::qe};
inline constexpr Species deuteron{"deuteron", 2.0 * constants::mp, constants::qe};
inline constexpr Species triton{"triton", 3.0 * constants::mp, constants::qe};
inline constexpr Species alpha{"alpha", 4.0 * constants::mp, 2.0 * constants::qe};

}  // namespace sim
