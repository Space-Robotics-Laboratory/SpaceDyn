cmake_minimum_required(VERSION 3.8)
project(spacedyn_ros)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Set package name
set(PACKAGE_NAME spacedyn_ros)

# Set ros packages dependencies
set(THIS_PACKAGE_INCLUDE_DEPENDS
  rclcpp
  ament_index_cpp
  urdf
  Eigen3
  geometry_msgs
  sensor_msgs
  std_msgs
  tf2
  tf2_ros
  control_msgs
  )

# Find dependencies
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

#### Set SpaceDyn folders
# Include directory
set(SPD_ROS_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Source directory
set(SPD_ROS_GEOMETRY_DIR 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/geometry)
set(SPD_ROS_LINKAGE_DIR 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/linkage)
set(SPD_ROS_MOTION_DIR 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/motion)
set(SPD_ROS_ROBOT_DIR 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/robot)
set(SPD_ROS_UTIL_DIR 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/util)

# Source files
set(SPD_ROS_GEOMETRY
  ${SPD_ROS_GEOMETRY_DIR}/accel.cpp
  ${SPD_ROS_GEOMETRY_DIR}/pose.cpp
  ${SPD_ROS_GEOMETRY_DIR}/transform.cpp
  ${SPD_ROS_GEOMETRY_DIR}/twist.cpp
  ${SPD_ROS_GEOMETRY_DIR}/wrench.cpp
  ${SPD_ROS_GEOMETRY_DIR}/inertia.cpp)
set(SPD_ROS_LINKAGE
  ${SPD_ROS_LINKAGE_DIR}/joint_state.cpp
  ${SPD_ROS_LINKAGE_DIR}/joint.cpp
  ${SPD_ROS_LINKAGE_DIR}/link_state.cpp
  ${SPD_ROS_LINKAGE_DIR}/link.cpp
  ${SPD_ROS_LINKAGE_DIR}/linkage.cpp)
set(SPD_ROS_MOTION
  ${SPD_ROS_MOTION_DIR}/kinematics.cpp
  ${SPD_ROS_MOTION_DIR}/integral.cpp
  ${SPD_ROS_MOTION_DIR}/dynamics.cpp)
set(SPD_ROS_ROBOT
  ${SPD_ROS_ROBOT_DIR}/model.cpp
  ${SPD_ROS_ROBOT_DIR}/robot.cpp
  ${SPD_ROS_ROBOT_DIR}/state_variable.cpp)
set(SPD_ROS_UTIL
  ${SPD_ROS_UTIL_DIR}/matrix_operation.cpp)


#### Make library to be included from other packages
# Add source
add_library(${PACKAGE_NAME} SHARED
  ${SPD_ROS_GEOMETRY}
  ${SPD_ROS_LINKAGE}
  ${SPD_ROS_MOTION}
  ${SPD_ROS_ROBOT}
  ${SPD_ROS_UTIL}
)

# Assign include directory
target_include_directories(${PACKAGE_NAME} PUBLIC
  $<BUILD_INTERFACE:${SPD_ROS_INCLUDE_DIRS}>
  $<INSTALL_INTERFACE:include>
)

# Install directories under ws/install/
install(DIRECTORY include/ DESTINATION include)
install(DIRECTORY
  model_sample
  DESTINATION share/${PROJECT_NAME})

# Link libraries used in spacedyn_ros
# And export dependencies
target_link_libraries(${PACKAGE_NAME} 
  ament_index_cpp::ament_index_cpp 
  Eigen3::Eigen)
ament_target_dependencies(${PACKAGE_NAME}
  ${THIS_PACKAGE_INCLUDE_DEPENDS}
)
ament_export_dependencies(${THIS_PACKAGE_INCLUDE_DEPENDS})
ament_export_targets(${PACKAGE_NAME} HAS_LIBRARY_TARGET)
ament_export_include_directories(${SPD_ROS_INCLUDE_DIRS})

install(
  TARGETS ${PACKAGE_NAME}
  EXPORT ${PACKAGE_NAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib
  INCLUDES DESTINATION include)
ament_package()
