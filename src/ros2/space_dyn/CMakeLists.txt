cmake_minimum_required(VERSION 3.8)
project(microgravity_manipulation)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(ament_index_cpp REQUIRED)
find_package(rclcpp REQUIRED)

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

##### START Add SpaceDyn #####
set(SPD_MAT_DIR spd/matrix/)
set(SPD_SPN_DIR spd/spn/)
set(SPD_SRC_DIR spd/src/)

set(SPD_MAT
  ${SPD_MAT_DIR}/matrix.cpp
  ${SPD_MAT_DIR}/vector.cpp)

set(SPD_SPN
  ${SPD_SPN_DIR}/rot.cpp
  ${SPD_SPN_DIR}/spn.cpp)

set(SPD_SRC
  ${SPD_SRC_DIR}/model.cpp
  ${SPD_SRC_DIR}/model_init.cpp
  ${SPD_SRC_DIR}/calc_SP.cpp
  ${SPD_SRC_DIR}/calc_SPN.cpp
  ${SPD_SRC_DIR}/calc_hh.cpp
  ${SPD_SRC_DIR}/calc_C.cpp
  ${SPD_SRC_DIR}/calc_Rg.cpp
  ${SPD_SRC_DIR}/calc_Je.cpp
  ${SPD_SRC_DIR}/calc_Jb.cpp
  ${SPD_SRC_DIR}/calc_GJe.cpp
  ${SPD_SRC_DIR}/calc_Lg.cpp
  ${SPD_SRC_DIR}/calc_JJ.cpp
  ${SPD_SRC_DIR}/calc_SS.cpp
  ${SPD_SRC_DIR}/calc_Xup_I.cpp
  ${SPD_SRC_DIR}/j_num.cpp
  ${SPD_SRC_DIR}/f_kin_e.cpp
  ${SPD_SRC_DIR}/f_kin_j.cpp
  ${SPD_SRC_DIR}/i_dyn.cpp
  ${SPD_SRC_DIR}/i_dyn_fix.cpp
  ${SPD_SRC_DIR}/f_dyn.cpp
  ${SPD_SRC_DIR}/f_dyn_fix.cpp
  ${SPD_SRC_DIR}/inner_force.cpp
  ${SPD_MAT_DIR}/matrix.cpp
  ${SPD_MAT_DIR}/vector.cpp
  ${SPD_SPN_DIR}/spn.cpp
  ${SPD_SPN_DIR}/rot.cpp
  ${SPD_SRC_DIR}/rref.cpp
  ${SPD_SRC_DIR}/aw.cpp
  ${SPD_SRC_DIR}/w2dQtn.cpp
  ${SPD_SRC_DIR}/f_dyn_rk.cpp
  ${SPD_SRC_DIR}/model_test.cpp
  ${SPD_SRC_DIR}/calc_GJb.cpp
  ${SPD_SRC_DIR}/calc_GJe_r.cpp)

##### END Add SpaceDyn #####
# Add Main 
set(SRC_DIR src)
set(SRC 
  ${SRC_DIR}/dual_arm_space_dyne.cpp)

# make library 
add_library(microgravity_manipulation SHARED 
              ${SPD_MAT} 
              ${SPD_SRC}
              ${SPD_SPN}
              ${SRC}) 

# Add dependencies of lib spd
target_link_libraries(microgravity_manipulation -lgsl -lgslcblas ament_index_cpp::ament_index_cpp)

ament_export_targets(microgravity_manipulation HAS_LIBRARY_TARGET)

# Add Test 
add_executable(test_gsl test/test_gsl.cpp)
add_executable(test_matrix_1 test/test_matrix_1.cpp)
add_executable(test_matrix_2 test/test_matrix_2.cpp)
add_executable(test_spd test/test_spd.cpp)
add_executable(test_dar test/test_dar.cpp)
target_link_libraries(test_gsl -lgsl -lgslcblas) # link gsl and cblass
target_link_libraries(test_matrix_1 microgravity_manipulation)
target_link_libraries(test_matrix_2 microgravity_manipulation)
target_link_libraries(test_spd microgravity_manipulation) 
target_link_libraries(test_dar microgravity_manipulation) 

install(TARGETS
  test_gsl
  test_matrix_1
  test_matrix_2
  test_spd
  test_dar
  DESTINATION lib/${PROJECT_NAME})

install(
  DIRECTORY include/${PROJECT_NAME}
  spd/
  DESTINATION include
)

install(
  TARGETS microgravity_manipulation
  EXPORT microgravity_manipulation
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION lib
  INCLUDES DESTINATION include
)

install(DIRECTORY
  model
  DESTINATION share/${PROJECT_NAME})

include_directories(
  include/${PROJECT_NAME}
  spd/
)

ament_package()
