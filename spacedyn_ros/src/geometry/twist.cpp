#include "spacedyn_ros/geometry/twist.hpp"
#include "eigen3/Eigen/Core"
#include "spacedyn_ros/geometry/pose.hpp"

#include "iostream"

namespace spacedyn_ros {

Twist::Twist(const Frame &frame, const Eigen::VectorXd &twist) {
  // Restrict twist size 6, as DOF in se3
  if (twist.size() != 6) {
    throw std::invalid_argument(
        "Error: Twist size is incorrect. Expect size is 6, Actual size is, " +
        std::to_string(twist.size()));
  }
  this->twist_ = twist;
  this->frame_ = frame;
}

Twist::Twist(const Frame &frame, const Eigen::Vector3d &linier_velocity,
             const Eigen::Vector3d &angular_velocity) {
  this->twist_.resize(6);
  this->twist_.block(0, 0, 3, 1) = linier_velocity;
  this->twist_.block(3, 0, 3, 1) = angular_velocity;
  this->frame_ = frame;
}

const Frame &Twist::getFrame() const { return frame_; }

const Eigen::VectorXd &Twist::getOriginTwist() const { return twist_; }
Eigen::Vector3d Twist::getOriginLinierVelocity() const { return twist_.block(0, 0, 3, 1); }
Eigen::Vector3d Twist::getOriginAngularVelocity() const { return twist_.block(3, 0, 3, 1); }

Twist Twist::operator+(const Twist &twist) const {
  if (frame_ != twist.getFrame()) {
    throw std::invalid_argument("Error: Twist must have the same frame");
  }
  return Twist(frame_, twist_ + twist.getOriginTwist());
}

Twist Twist::operator-(const Twist &twist) const {
  if (frame_ != twist.getFrame()) {
    throw std::invalid_argument("Error: Twist must have the same frame");
  }
  return Twist(frame_, twist_ - twist.getOriginTwist());
}

Twist Twist::computePointTwist(const Transform &tf_to_point) const {
  Eigen::VectorXd twist_point_in_world(6);
  Eigen::Vector3d linier = getOriginLinierVelocity();
  Eigen::Vector3d angular = getOriginAngularVelocity();

  if (frame_ != tf_to_point.getFrame()) {
    throw std::invalid_argument("Error: Twist and Transform should have the same Frame");
  }

  if (tf_to_point.getFrame() == Frame::kLocal) {
    twist_point_in_world.block(0, 0, 3, 1) =
        linier + angular.cross(tf_to_point.getTransform().translation());
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
  Eigen::VectorXd twist_point_in_world(6);
  Eigen::Vector3d linier = getOriginLinierVelocity();
  Eigen::Vector3d angular = getOriginAngularVelocity();
  twist_point_in_world.block(0, 0, 3, 1) =
      linier + angular.cross(translation_to_point_in_world_frame);
  twist_point_in_world.block(3, 0, 3, 1) = angular;

  return Twist(Frame::kWorld, twist_point_in_world);
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
