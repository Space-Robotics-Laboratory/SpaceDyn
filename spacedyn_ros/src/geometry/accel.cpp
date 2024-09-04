#include "spacedyn_ros/geometry/accel.hpp"
#include "eigen3/Eigen/Core"

#include "iostream"

namespace spacedyn_ros {

Accel::Accel(const Frame &frame, const Eigen::VectorXd &accel) {
  // Restrict accel size 6, as DOF in se3
  if (accel.size() != 6) {
    throw std::invalid_argument(
        "Error: Accel size is incorrect. Expect size is 6, Actual size is, " +
        std::to_string(accel.size()));
  }
  this->frame_ = frame;
  this->accel_ = accel;
}

const Frame &Accel::getFrame() const { return frame_; }
const Eigen::VectorXd &Accel::getOriginAccel() const { return accel_; }
Eigen::Vector3d Accel::getOriginLinierAcceleration() const { return accel_.block(0, 0, 3, 1); }
Eigen::Vector3d Accel::getOriginAngularAcceleration() const { return accel_.block(3, 0, 3, 1); }

Accel Accel::operator+(const Accel &accel) const {
  if (frame_ != accel.getFrame()) {
    throw std::invalid_argument("Error: Accel must have the same frame");
  }
  return Accel(frame_, accel_ + accel.getOriginAccel());
}

Accel Accel::operator-(const Accel &accel) const {
  if (frame_ != accel.getFrame()) {
    throw std::invalid_argument("Error: Accel must have the same frame");
  }
  return Accel(frame_, accel_ - accel.getOriginAccel());
}

Accel Accel::computePointAccel(const Twist &twist,
                               const Eigen::Vector3d &translation_to_point_in_world_frame) const {
  if (frame_ != Frame::kWorld) {
    throw std::invalid_argument("Error: Accel must be in World frame");
  }
  if (twist.getFrame() != Frame::kWorld) {
    throw std::invalid_argument("Error: Twist must be in World frame");
  }

  Eigen::VectorXd accel_point_in_space(6);
  Eigen::Vector3d linier = getOriginLinierAcceleration();
  Eigen::Vector3d angular = getOriginAngularAcceleration();
  auto omega = twist.getOriginAngularVelocity();

  accel_point_in_space.block(0, 0, 3, 1) =
      linier + angular.cross(translation_to_point_in_world_frame) +
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
