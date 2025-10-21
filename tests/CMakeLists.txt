find_package(nlohmann_json CONFIG QUIET)
if(NOT nlohmann_json_FOUND)
  FetchContent_Declare(
    nlohmann_json
    URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
  )
  FetchContent_MakeAvailable(nlohmann_json)
endif()

add_executable(
  tokamak_sim_tests
  tests_smoke.cpp
  tests_regression.cpp
)

target_link_libraries(
  tokamak_sim_tests
  PRIVATE
    tokamak_sim_core
    gtest_main
    nlohmann_json::nlohmann_json
)

target_compile_definitions(
  tokamak_sim_tests
  PRIVATE
    TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data"
)

add_test(NAME tokamak_sim_tests COMMAND tokamak_sim_tests)
