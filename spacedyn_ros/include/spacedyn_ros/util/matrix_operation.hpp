#ifndef SPACEDYN_ROS_GEOMETRY_UTIL_HPP_
#define SPACEDYN_ROS_GEOMETRY_UTIL_HPP_
#include "eigen3/Eigen/Core"

namespace spacedyn_ros {

Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d &v);
bool isPositiveDefinite(const Eigen::MatrixXd &m);
bool isSymmetric(const Eigen::MatrixXd &m);

} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_GEOMETRY_UTIL_HPP_
