cmake_minimum_required(VERSION 3.8)
project(spd_example)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Set package name
set(PACKAGE_NAME spd_example)

# Set ros packages dependencies
set(THIS_PACKAGE_INCLUDE_DEPENDS
  spacedyn_ros
)

# find dependencies
find_package(ament_cmake REQUIRED)
foreach(Dependency IN ITEMS ${THIS_PACKAGE_INCLUDE_DEPENDS})
  find_package(${Dependency} REQUIRED)
endforeach()

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

add_executable(zhanibekov_effect_sim src/zhanibekov_effect_sim.cpp)
target_link_libraries(zhanibekov_effect_sim spacedyn_ros::spacedyn_ros)
ament_target_dependencies(zhanibekov_effect_sim spacedyn_ros)

target_include_directories(zhanibekov_effect_sim
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Install Shared Directories
set(INSTALL_SHARE_DIRS
  launch
  config
  urdf
)
foreach(DIR IN LISTS INSTALL_SHARE_DIRS)
  install(DIRECTORY ${DIR} DESTINATION share/${PROJECT_NAME})
endforeach()

install(TARGETS
  zhanibekov_effect_sim
  DESTINATION lib/${PROJECT_NAME}
)

ament_package()
