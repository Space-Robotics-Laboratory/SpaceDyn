#ifndef SPACEDYN_ROS_GEOMETRY_UTIL_HPP_
#define SPACEDYN_ROS_GEOMETRY_UTIL_HPP_
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

namespace Eigen {
typedef Eigen::Matrix<double, 6, 1> Vector6d;
}

namespace spacedyn_ros {

Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d &v);
Eigen::Quaterniond quatFromTwoVectors(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2);
bool isPositiveDefinite(const Eigen::MatrixXd &m);
bool isSymmetric(const Eigen::MatrixXd &m);

} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_GEOMETRY_UTIL_HPP_
