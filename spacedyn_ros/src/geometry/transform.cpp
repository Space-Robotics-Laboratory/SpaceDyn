#include "spacedyn_ros/geometry/transform.hpp"
#include "spacedyn_ros/geometry/frame.hpp"
#include <iostream>

namespace spacedyn_ros {
Transform::Transform(const Frame frame, const Eigen::Isometry3d &transform) {
  // T = [R, p; 0, 1] = translation x rotation
  try {
    checkRotationMatrixNormalized(transform.rotation());
    this->frame_ = frame;
    this->transform_ = transform;
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw std::invalid_argument("Error: Failed to create transform. ");
  }
}

Transform::Transform(const Frame frame, const Eigen::Matrix3d &rotation,
                     const Eigen::Vector3d &translation) {
  // T = [R, p; 0, 1] = translation x rotation
  try {
    checkRotationMatrixNormalized(rotation);
    this->frame_ = frame;
    this->transform_ = Eigen::Translation3d(translation) * Eigen::Quaterniond(rotation);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw std::invalid_argument("Error: Failed to create transform. ");
  }
}

Transform::Transform(const Frame frame, const Eigen::Quaterniond &rotation,
                     const Eigen::Vector3d &translation) {
  // T = [R, p; 0, 1] = translation x rotation
  try {
    checkQuaternionNormalized(rotation);
    this->frame_ = frame;
    this->transform_ = Eigen::Translation3d(translation) * rotation;
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw std::invalid_argument("Error: Failed to create transform. ");
  }
}

void Transform::checkRotationMatrixNormalized(const Eigen::Matrix3d &rotation_matrix) const {
  if (std::abs(rotation_matrix.determinant() - 1.0) > 1e-6) {
    throw std::runtime_error("Error: Rotation matrix is not normalized.");
  }
}

void Transform::checkQuaternionNormalized(const Eigen::Quaterniond &quaternion) const {
  if (std::abs(quaternion.norm() - 1.0) > 1e-6) {
    throw std::runtime_error("Error: Quaternion is not normalized.");
  }
}

const Frame &Transform::getFrame() const { return frame_; }
const Eigen::Isometry3d &Transform::getTransform() const { return transform_; }
Transform Transform::inverse() const {
  return Transform(frame_, transform_.inverse().rotation(), transform_.inverse().translation());
}

Transform Transform::transform(const Transform &affected_tf) const {
  if (frame_ != affected_tf.getFrame()) {
    throw std::runtime_error("Error: Frame mismatch. Cannot transform between different frames.");
  }
  if (frame_ == Frame::kWorld) {
    // T' = Ts * T
    return Transform(frame_, transform_ * affected_tf.getTransform());
  } else {
    // T' = T * Tb
    return Transform(frame_, affected_tf.getTransform() * transform_);
  }
}

Eigen::Vector3d Transform::getTranslation() const { return transform_.translation(); }
Eigen::Matrix3d Transform::getRotation() const { return transform_.rotation(); }

} // namespace spacedyn_ros
