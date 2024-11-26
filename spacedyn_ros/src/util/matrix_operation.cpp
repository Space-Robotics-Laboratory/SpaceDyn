#include "spacedyn_ros/util/matrix_operation.hpp"
#include <eigen3/Eigen/Cholesky>

namespace spacedyn_ros {
Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d &v) {
  Eigen::Matrix3d skew_symmetric;
  skew_symmetric << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
  return skew_symmetric;
}

Eigen::Quaterniond quatFromTwoVectors(const Eigen::Vector3d &a, const Eigen::Vector3d &b) {
  // Normalize both vectors
  Eigen::Vector3d a_norm = a.normalized();
  Eigen::Vector3d b_norm = b.normalized();

  // Compute the axis of rotation (cross product)
  Eigen::Vector3d axis = a_norm.cross(b_norm);

  // Compute the angle between the vectors (dot product and arccosine)
  double cos_theta = a_norm.dot(b_norm);
  double angle = std::acos(cos_theta);

  // If the vectors are already aligned, return the identity quaternion
  if (cos_theta > 1.0 - 1e-6) {
    return Eigen::Quaterniond::Identity();
  }

  // If the vectors are opposite, find an orthogonal vector to define the rotation axis
  if (cos_theta < -1.0 + 1e-6) {
    Eigen::Vector3d orthogonal = Eigen::Vector3d::UnitX().cross(a_norm);
    if (orthogonal.norm() < 1e-6) {
      orthogonal = Eigen::Vector3d::UnitY().cross(a_norm);
    }
    orthogonal.normalize();
    return Eigen::Quaterniond(Eigen::AngleAxisd(M_PI, orthogonal));
  }

  // Normalize the axis of rotation
  axis.normalize();

  // Construct the quaternion from the axis and angle
  return Eigen::Quaterniond(Eigen::AngleAxisd(angle, axis));
}

bool isPositiveDefinite(const Eigen::MatrixXd &mat) {
  Eigen::LLT<Eigen::MatrixXd> llt(mat); // Cholesky decomposition
  return llt.info() == Eigen::Success;  // Successful decomposition means positive definite
}

bool isSymmetric(const Eigen::MatrixXd &mat) { return mat.isApprox(mat.transpose()); }
} // namespace spacedyn_ros
