#include "spacedyn_ros/util/matrix_operation.hpp"
#include "eigen3/Eigen/Cholesky"

namespace spacedyn_ros {
Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d &v) {
  Eigen::Matrix3d skew_symmetric;
  skew_symmetric << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
  return skew_symmetric;
}

bool isPositiveDefinite(const Eigen::MatrixXd &mat) {
  Eigen::LLT<Eigen::MatrixXd> llt(mat); // Cholesky decomposition
  return llt.info() == Eigen::Success;  // Successful decomposition means positive definite
}

bool isSymmetric(const Eigen::MatrixXd &mat) { return mat.isApprox(mat.transpose()); }
} // namespace spacedyn_ros
