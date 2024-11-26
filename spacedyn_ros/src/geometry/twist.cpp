#include "spacedyn_ros/geometry/twist.hpp"
#include "spacedyn_ros/geometry/pose.hpp"
#include "spacedyn_ros/util/matrix_operation.hpp"
#include <eigen3/Eigen/Core>
#include <iostream>

namespace spacedyn_ros {

Twist::Twist(const Frame &frame, const Eigen::Vector6d &twist) {
  // Restrict twist size 6, as DOF in se3
  if (twist.size() != 6) {
    throw std::invalid_argument(
        "Error: Twist size is incorrect. Expect size is 6, Actual size is, " +
        std::to_string(twist.size()));
  }
  this->twist_ = twist;
  this->frame_ = frame;
}

Twist::Twist(const Frame &frame, const Eigen::Vector3d &linear_velocity,
             const Eigen::Vector3d &angular_velocity) {
  this->twist_.resize(6);
  this->twist_.block(0, 0, 3, 1) = linear_velocity;
  this->twist_.block(3, 0, 3, 1) = angular_velocity;
  this->frame_ = frame;
}

const Frame &Twist::getFrame() const { return frame_; }

const Eigen::Vector6d &Twist::getTwist() const { return twist_; }
Eigen::Vector3d Twist::getLinearVelocity() const { return twist_.block(0, 0, 3, 1); }
Eigen::Vector3d Twist::getAngularVelocity() const { return twist_.block(3, 0, 3, 1); }

Twist Twist::operator+(const Twist &twist) const {
  if (frame_ != twist.getFrame()) {
    throw std::invalid_argument("Error: Twist must have the same frame");
  }
  return Twist(frame_, twist_ + twist.getTwist());
}

Twist Twist::operator-(const Twist &twist) const {
  if (frame_ != twist.getFrame()) {
    throw std::invalid_argument("Error: Twist must have the same frame");
  }
  return Twist(frame_, twist_ - twist.getTwist());
}

Twist Twist::getTwistInFrame(const Frame &frame, const Pose &pose) const {
  if (frame == frame_) {
    // No need to change the frame
    return Twist(frame_, twist_);
  } else if (frame == Frame::kLocal) { // frame_ == Frame::kWorld
    // Change the frame to local
    // v_a = Ras * v_s
    Eigen::Vector3d linear_velocity =
        pose.getAttitudeInWorldFrame().transpose() * getLinearVelocity();
    Eigen::Vector3d angular_velocity =
        pose.getAttitudeInWorldFrame().transpose() * getAngularVelocity();
    return Twist(Frame::kLocal, linear_velocity, angular_velocity);
  } else if (frame == Frame::kWorld) { // frame_ == Frame::kLocal
    // Change the frame to world
    // v_s = Rsa * v_a
    Eigen::Vector3d linear_velocity = pose.getAttitudeInWorldFrame() * getLinearVelocity();
    Eigen::Vector3d angular_velocity = pose.getAttitudeInWorldFrame() * getAngularVelocity();
    return Twist(Frame::kWorld, linear_velocity, angular_velocity);
  }
  throw std::invalid_argument("Error: Unknown frame type. ");
}

Twist Twist::computePointTwist(const Transform &tf_to_point) const {
  Eigen::Vector6d twist_point_in_world;
  Eigen::Vector3d linear = getLinearVelocity();
  Eigen::Vector3d angular = getAngularVelocity();

  if (frame_ != tf_to_point.getFrame()) {
    throw std::invalid_argument("Error: Twist and Transform should have the same Frame");
  }

  if (tf_to_point.getFrame() == Frame::kLocal) {
    twist_point_in_world.block(0, 0, 3, 1) =
        linear + angular.cross(tf_to_point.getTransform().translation());
    twist_point_in_world.block(3, 0, 3, 1) = angular;
  } else if (tf_to_point.getFrame() == Frame::kWorld) {
    // TODO: implement world
    throw std::runtime_error("Error: computePointTwist is not implemented for world frame");
  }
  return Twist(Frame::kLocal, twist_point_in_world);
}

Twist Twist::computePointTwist(const Eigen::Vector3d &translation_to_point_in_world_frame) const {
  if (frame_ != Frame::kWorld) {
    throw std::invalid_argument("Error: Twist must be in World frame to use this function");
  }
  Eigen::Vector6d twist_point_in_world;
  Eigen::Vector3d linear = getLinearVelocity();
  Eigen::Vector3d angular = getAngularVelocity();
  twist_point_in_world.block(0, 0, 3, 1) =
      linear + angular.cross(translation_to_point_in_world_frame);
  twist_point_in_world.block(3, 0, 3, 1) = angular;

  return Twist(Frame::kWorld, twist_point_in_world);
}

Eigen::Quaterniond Twist::computeDerivativeAttitude(const Pose &pose) const {
  Eigen::Quaterniond q = pose.getQuaternionInWorldFrame();
  Eigen::Vector3d omega_local;
  Eigen::Quaterniond q_dot;
  switch (frame_) {
  case Frame::kWorld:
    omega_local = pose.getAttitudeInWorldFrame().inverse() * getAngularVelocity();
    break;

  case Frame::kLocal:
    // q_dot = 0.5 * q * [0, omega]
    omega_local = getAngularVelocity();
    break;

  default:
    throw std::invalid_argument("Error: Unknown frame type. ");
  }

  q_dot =
      q * Eigen::Quaterniond(0, 0.5 * omega_local[0], 0.5 * omega_local[1], 0.5 * omega_local[2]);

  return q_dot;
}

geometry_msgs::msg::Twist Twist::toRosMessage() const {
  geometry_msgs::msg::Twist msg;
  msg.linear.x = twist_(0);
  msg.linear.y = twist_(1);
  msg.linear.z = twist_(2);
  msg.angular.x = twist_(3);
  msg.angular.y = twist_(4);
  msg.angular.z = twist_(5);
  return msg;
}
} // namespace spacedyn_ros
