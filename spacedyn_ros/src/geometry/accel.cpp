#include "spacedyn_ros/geometry/accel.hpp"
#include "spacedyn_ros/util/matrix_operation.hpp"
#include <eigen3/Eigen/Core>
#include <iostream>

namespace spacedyn_ros {

Accel::Accel(const Frame &frame, const Eigen::Vector6d &accel) {
  // Restrict accel size 6, as DOF in se3
  if (accel.size() != 6) {
    throw std::invalid_argument(
        "Error: Accel size is incorrect. Expect size is 6, Actual size is, " +
        std::to_string(accel.size()));
  }
  this->frame_ = frame;
  this->accel_ = accel;
}

Accel::Accel(const Frame &frame, const Eigen::Vector3d &linear_acceleration,
             const Eigen::Vector3d &angular_acceleration) {
  this->frame_ = frame;
  this->accel_.block(0, 0, 3, 1) = linear_acceleration;
  this->accel_.block(3, 0, 3, 1) = angular_acceleration;
}

const Frame &Accel::getFrame() const { return frame_; }
const Eigen::Vector6d &Accel::getAccel() const { return accel_; }
Eigen::Vector3d Accel::getLinearAcceleration() const { return accel_.block(0, 0, 3, 1); }
Eigen::Vector3d Accel::getAngularAcceleration() const { return accel_.block(3, 0, 3, 1); }

Accel Accel::operator+(const Accel &accel) const {
  if (frame_ != accel.getFrame()) {
    throw std::invalid_argument("Error: Accel must have the same frame");
  }
  return Accel(frame_, accel_ + accel.getAccel());
}

Accel Accel::operator-(const Accel &accel) const {
  if (frame_ != accel.getFrame()) {
    throw std::invalid_argument("Error: Accel must have the same frame");
  }
  return Accel(frame_, accel_ - accel.getAccel());
}

Accel Accel::getAccelInFrame(const Frame &frame, const Pose &pose) const {
  if (frame == frame_) {
    // No need to change the frame
    return Accel(frame_, accel_);
  } else if (frame == Frame::kLocal) { // frame_ == Frame::kWorld
    // Change the frame to local
    // a_a = Ras * a_s
    Eigen::Vector3d linear_acceleration =
        pose.getAttitudeInWorldFrame().transpose() * getLinearAcceleration();
    Eigen::Vector3d angular_acceleration =
        pose.getAttitudeInWorldFrame().transpose() * getAngularAcceleration();
    return Accel(Frame::kLocal, linear_acceleration, angular_acceleration);
  } else if (frame == Frame::kWorld) { // frame_ == Frame::kLocal
    // Change the frame to world
    // a_s = Rsa * a_a
    Eigen::Vector3d linear_acceleration = pose.getAttitudeInWorldFrame() * getLinearAcceleration();
    Eigen::Vector3d angular_acceleration =
        pose.getAttitudeInWorldFrame() * getAngularAcceleration();
    return Accel(Frame::kWorld, linear_acceleration, angular_acceleration);
  }
  throw std::invalid_argument("Error: Unknown frame type. ");
}

Accel Accel::computePointAccel(const Twist &twist,
                               const Eigen::Vector3d &translation_to_point_in_world_frame) const {
  if (frame_ != Frame::kWorld) {
    throw std::invalid_argument("Error: Accel must be in World frame");
  }
  if (twist.getFrame() != Frame::kWorld) {
    throw std::invalid_argument("Error: Twist must be in World frame");
  }

  Eigen::Vector6d accel_point_in_space(6);
  Eigen::Vector3d linear = getLinearAcceleration();
  Eigen::Vector3d angular = getAngularAcceleration();
  auto omega = twist.getAngularVelocity();

  accel_point_in_space.block(0, 0, 3, 1) =
      linear + angular.cross(translation_to_point_in_world_frame) +
      omega.cross(omega.cross(translation_to_point_in_world_frame));
  accel_point_in_space.block(3, 0, 3, 1) = angular;
  return Accel(frame_, accel_point_in_space);
}

geometry_msgs::msg::Accel Accel::toRosMessage() const {
  geometry_msgs::msg::Accel msg;
  msg.linear.x = accel_(0);
  msg.linear.y = accel_(1);
  msg.linear.z = accel_(2);
  msg.angular.x = accel_(3);
  msg.angular.y = accel_(4);
  msg.angular.z = accel_(5);
  return msg;
}

} // namespace spacedyn_ros
