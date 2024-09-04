#include "spacedyn_ros/geometry/pose.hpp"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"
#include "iostream"
#include "spacedyn_ros/geometry/twist.hpp"

namespace spacedyn_ros {

Pose::Pose(const Eigen::Isometry3d &pose) {
  try {
    checkRotationMatrixNormalized(pose.rotation());
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw std::invalid_argument("Error: Failed to create pose. ");
  }
  pose_in_world_ = pose;
}
Pose::Pose(const Eigen::Vector3d &position, const Eigen::Matrix3d &attitude) {
  // The order of the following Trans() and Rotation() is important.
  try {
    checkRotationMatrixNormalized(attitude);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw std::invalid_argument("Error: Failed to create pose. ");
  }
  pose_in_world_ = Eigen::Translation3d(position) * Eigen::Quaterniond(attitude);
}
Pose::Pose(const Eigen::Vector3d &position, const Eigen::Quaterniond &attitude) {
  try {
    checkQuaternionNormalized(attitude);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw std::invalid_argument("Error: Failed to create pose. ");
  }
  pose_in_world_ = Eigen::Translation3d(position) * attitude;
}

void Pose::checkRotationMatrixNormalized(const Eigen::Matrix3d &rotation_matrix) const {
  if (std::abs(rotation_matrix.determinant() - 1.0) > 1e-6) {
    throw std::invalid_argument("Error: Rotation matrix violate norm. ");
  }
}

void Pose::checkQuaternionNormalized(const Eigen::Quaterniond &quaternion) const {
  if (std::abs(quaternion.norm() - 1.0) > 1e-6) {
    throw std::invalid_argument("Error: Quaternion violate norm. ");
  }
}

const Eigen::Isometry3d &Pose::getOriginPose() const { return pose_in_world_; }
Eigen::Vector3d Pose::getOriginPosition() const { return pose_in_world_.translation(); }
Eigen::MatrixXd Pose::getOriginAttitude() const { return pose_in_world_.rotation(); }
Eigen::Quaterniond Pose::getOriginQuaternion() const {
  return Eigen::Quaterniond(pose_in_world_.rotation());
}

Eigen::Vector3d Pose::computeTranslationToPoint(const Pose &point_pose) const {
  Eigen::Vector3d translation_to_point = point_pose.getOriginPosition() - getOriginPosition();
  return translation_to_point;
}

Transform Pose::computeTransformToPoint(const Frame &frame, const Pose &point_pose) const {
  Transform tf_to_point;
  Eigen::Isometry3d ism_to_point;
  switch (frame) {
  case Frame::kWorld:
    ism_to_point = point_pose.getOriginPose() * pose_in_world_.inverse();
    tf_to_point = Transform(Frame::kWorld, ism_to_point.rotation(), ism_to_point.translation());
    break;

  case Frame::kLocal:
    ism_to_point = pose_in_world_.inverse() * point_pose.getOriginPose();
    tf_to_point = Transform(Frame::kLocal, ism_to_point.rotation(), ism_to_point.translation());
    break;
  }

  return tf_to_point;
}

Pose Pose::computePointPose(const Transform &tf_to_point) const {
  auto frame = tf_to_point.getFrame();
  Eigen::Isometry3d pose_at_point;
  switch (frame) {
  case Frame::kWorld:
    // P' = T_w * P
    pose_at_point = tf_to_point.getTransform() * pose_in_world_;
    break;

  case Frame::kLocal:
    // P' = P * T_w
    pose_at_point = pose_in_world_ * tf_to_point.getTransform();
    break;
  }
  return Pose(pose_at_point);
}

Pose Pose::computeOriginPoseFromPointPose(const Transform &tf_to_point,
                                          const Pose &pose_at_point) const {
  auto frame = tf_to_point.getFrame();
  Eigen::Isometry3d pose_at_origin;
  switch (frame) {
  case Frame::kWorld:
    // P = T_w^-1 * P'
    pose_at_origin = tf_to_point.getTransform().inverse() * pose_at_point.getOriginPose();
    break;

  case Frame::kLocal:
    // P = P' * T_b^-1
    pose_at_origin = pose_at_point.getOriginPose() * tf_to_point.getTransform().inverse();
    break;
  }
  return Pose(pose_at_origin);
}

Twist Pose::computeTwistInLocalFrame(const Twist &twist_in_world_frame) const {
  Eigen::Vector3d linier_velocity;
  Eigen::Vector3d angular_velocity;
  switch (twist_in_world_frame.getFrame()) {
  case Frame::kWorld:
    linier_velocity =
        pose_in_world_.rotation().inverse() * twist_in_world_frame.getOriginLinierVelocity();
    angular_velocity =
        pose_in_world_.rotation().inverse() * twist_in_world_frame.getOriginAngularVelocity();
    break;

  case Frame::kLocal:
    linier_velocity = twist_in_world_frame.getOriginLinierVelocity();
    angular_velocity = twist_in_world_frame.getOriginAngularVelocity();
    break;

  default:
    throw std::invalid_argument("Error: Unknown frame type. ");
  }

  return Twist(Frame::kLocal, linier_velocity, angular_velocity);
}

Wrench Pose::computeWrenchInWorldFrame(const Wrench &wrench_in_local_frame) const {
  if (wrench_in_local_frame.getFrame() != Frame::kLocal) {
    throw std::invalid_argument("Error: wrench should be in local frame. ");
  }
  // Ws = Rsa * Wa
  Eigen::Vector3d force = pose_in_world_.rotation() * wrench_in_local_frame.getOriginForce();
  Eigen::Vector3d torque = pose_in_world_.rotation() * wrench_in_local_frame.getOriginTorque();
  return Wrench(Frame::kWorld, force, torque);
}

Eigen::Quaterniond Pose::computeDerivativeAttitude(const Twist &twist) const {
  Eigen::Quaterniond q = Eigen::Quaterniond(pose_in_world_.rotation());
  Eigen::Vector3d omega_local;
  Eigen::Quaterniond q_dot;
  switch (twist.getFrame()) {
  case Frame::kWorld:
    omega_local = pose_in_world_.rotation().inverse() * twist.getOriginAngularVelocity();
    break;

  case Frame::kLocal:
    // q_dot = 0.5 * q * [0, omega]
    omega_local = twist.getOriginAngularVelocity();
    break;

  default:
    throw std::invalid_argument("Error: Unknown frame type. ");
  }

  q_dot =
      q * Eigen::Quaterniond(0, 0.5 * omega_local[0], 0.5 * omega_local[1], 0.5 * omega_local[2]);

  return q_dot;
}

geometry_msgs::msg::Pose Pose::toRosMessage() const {
  geometry_msgs::msg::Pose pose_msg;
  Eigen::Vector3d position = getOriginPosition();
  Eigen::Quaterniond attitude(pose_in_world_.rotation());
  pose_msg.position.x = position.x();
  pose_msg.position.y = position.y();
  pose_msg.position.z = position.z();
  pose_msg.orientation.x = attitude.x();
  pose_msg.orientation.y = attitude.y();
  pose_msg.orientation.z = attitude.z();
  pose_msg.orientation.w = attitude.w();
  return pose_msg;
}

geometry_msgs::msg::TransformStamped Pose::toRosMessage(const std::string &frame_name,
                                                        const std::string &child_frame_name) const {
  geometry_msgs::msg::TransformStamped msg;
  msg.header.frame_id = frame_name;
  msg.child_frame_id = child_frame_name;
  auto position = getOriginPosition();
  auto attitude = getOriginQuaternion();
  msg.transform.translation.x = position.x();
  msg.transform.translation.y = position.y();
  msg.transform.translation.z = position.z();
  msg.transform.rotation.x = attitude.x();
  msg.transform.rotation.y = attitude.y();
  msg.transform.rotation.z = attitude.z();
  msg.transform.rotation.w = attitude.w();
  return msg;
}
} // namespace spacedyn_ros
