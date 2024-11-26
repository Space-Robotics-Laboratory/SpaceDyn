#ifndef SPACEDYN_ROS_TRANSFORM_HPP_
#define SPACEDYN_ROS_TRANSFORM_HPP_

#include "spacedyn_ros/geometry/frame.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <geometry_msgs/msg/transform.hpp>

namespace spacedyn_ros {
class Transform {
private:
  Eigen::Isometry3d transform_;
  Frame frame_;
  void checkRotationMatrixNormalized(const Eigen::Matrix3d &rotation_matrix) const;
  void checkQuaternionNormalized(const Eigen::Quaterniond &quaternion) const;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Transform(const Frame frame = Frame::kLocal,
            const Eigen::Isometry3d &transform = Eigen::Isometry3d::Identity());
  Transform(const Frame frame, const Eigen::Matrix3d &rotation, const Eigen::Vector3d &translation);
  Transform(const Frame frame, const Eigen::Quaterniond &rotation,
            const Eigen::Vector3d &translation);
  ~Transform() = default;

  Transform inverse() const;
  Transform transform(const Transform &affected_tf) const;

  const Frame &getFrame() const;
  const Eigen::Isometry3d &getTransform() const;
  Eigen::Vector3d getTranslation() const;
  Eigen::Matrix3d getRotation() const;
};

} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_TRANSFORM_HPP_
