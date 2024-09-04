#include "spacedyn_ros/geometry/wrench.hpp"
#include "eigen3/Eigen/Core"

#include "iostream"

namespace spacedyn_ros {

Wrench::Wrench(const Frame &frame, const Eigen::VectorXd &wrench) {
  // Restrict wrench size 6, as DOF in se3
  if (wrench.size() != 6) {
    throw std::invalid_argument(
        "Error: Wrench size is incorrect. Expect size is 6, Actual size is, " +
        std::to_string(wrench.size()));
  }
  this->wrench_ = wrench;
  this->frame_ = frame;
}

Wrench::Wrench(const Frame &frame, const Eigen::Vector3d &force, const Eigen::Vector3d &torque) {
  this->frame_ = frame;
  this->wrench_.resize(6);
  this->wrench_.block(0, 0, 3, 1) = force;
  this->wrench_.block(3, 0, 3, 1) = torque;
}

const Frame &Wrench::getFrame() const { return frame_; }

const Eigen::VectorXd &Wrench::getOriginWrench() const { return wrench_; }

Eigen::Vector3d Wrench::getOriginForce() const { return wrench_.block(0, 0, 3, 1); }

Eigen::Vector3d Wrench::getOriginTorque() const { return wrench_.block(3, 0, 3, 1); }

Wrench Wrench::operator+(const Wrench &wrench) const {
  if (frame_ != wrench.getFrame()) {
    throw std::invalid_argument("Error: Wrench must have the same frame");
  }
  return Wrench(frame_, wrench_ + wrench.getOriginWrench());
}

Wrench Wrench::operator-() const { return Wrench(frame_, -wrench_); }

Wrench Wrench::operator-(const Wrench &wrench) const {
  if (frame_ != wrench.getFrame()) {
    throw std::invalid_argument("Error: Wrench must have the same frame");
  }
  return Wrench(frame_, wrench_ - wrench.getOriginWrench());
}

Wrench Wrench::operator+=(const Wrench &wrench) {
  if (frame_ != wrench.getFrame()) {
    throw std::invalid_argument("Error: Wrench must have the same frame");
  }
  wrench_ += wrench.getOriginWrench();
  return *this;
}

Wrench Wrench::computePointWrench(const Transform &tf_to_point) const {
  Eigen::VectorXd wrench_point_in_world(6);
  Eigen::Vector3d force = getOriginForce();
  Eigen::Vector3d torque = getOriginTorque();

  if (frame_ != tf_to_point.getFrame()) {
    throw std::invalid_argument("Error: Wrench and Transform should have the same Frame");
  }

  if (tf_to_point.getFrame() == Frame::kLocal) {
    wrench_point_in_world.block(0, 0, 3, 1) = force;
    wrench_point_in_world.block(3, 0, 3, 1) =
        torque + force.cross(tf_to_point.getTransform().translation());
  } else if (tf_to_point.getFrame() == Frame::kWorld) {
    // TODO: implement world
    throw std::runtime_error("Error: computePointWrench is not implemented for world frame");
  }
  return Wrench(Frame::kLocal, wrench_point_in_world);
}

Wrench Wrench::computeOriginWrenchFromPointWrench(
    const Eigen::Vector3d &translation_to_point_in_world_frame) const {
  if (frame_ != Frame::kWorld) {
    throw std::invalid_argument("Error: Wrench must be in World frame to use this function");
  }
  Eigen::VectorXd wrench_point_in_world(6);
  Eigen::Vector3d force = getOriginForce();
  Eigen::Vector3d torque = getOriginTorque();
  wrench_point_in_world.block(0, 0, 3, 1) = force;
  wrench_point_in_world.block(3, 0, 3, 1) =
      torque + translation_to_point_in_world_frame.cross(force);
  return Wrench(Frame::kWorld, wrench_point_in_world);
}

geometry_msgs::msg::Wrench Wrench::toRosMessage() const {
  geometry_msgs::msg::Wrench msg;
  msg.force.x = wrench_(0);
  msg.force.y = wrench_(1);
  msg.force.z = wrench_(2);
  msg.torque.x = wrench_(3);
  msg.torque.y = wrench_(4);
  msg.torque.z = wrench_(5);
  return msg;
}

} // namespace spacedyn_ros
