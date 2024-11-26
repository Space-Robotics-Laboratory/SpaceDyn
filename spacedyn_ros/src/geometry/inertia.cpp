#include "spacedyn_ros/geometry/inertia.hpp"
#include "spacedyn_ros/util/matrix_operation.hpp"
#include <iostream>

namespace spacedyn_ros {
Inertia::Inertia(const Frame &frame, const double mass, const Eigen::Matrix3d &inertia) {
  try {
    checkMass(mass);
    checkInertia(inertia);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Inertia initialization failed. ");
  }

  this->frame_ = frame;
  this->mass_ = mass;
  this->inertia_ = inertia;
}

Inertia Inertia::getInertiaInFrame(const Frame &frame, const Pose &pose) const {
  Eigen::Matrix3d inertia_in_frame;
  if (frame == frame_) {
    // No need to change the frame
    return Inertia(frame_, mass_, inertia_);
  } else if (frame == Frame::kLocal) { // frame_ == Frame::kWorld
    // Change the frame to local
    // Ia = Rsa^T * Is * Rsa
    Eigen::Matrix3d inertia_in_local =
        pose.getAttitudeInWorldFrame().transpose() * inertia_ * pose.getAttitudeInWorldFrame();
    return Inertia(Frame::kLocal, mass_, inertia_in_local);
  } else if (frame == Frame::kWorld) { // frame_ == Frame::kLocal
    // Change the frame to world
    // Is = Rsa * Ia * Rsa^T
    Eigen::Matrix3d inertia_in_world =
        pose.getAttitudeInWorldFrame() * inertia_ * pose.getAttitudeInWorldFrame().transpose();
    return Inertia(Frame::kWorld, mass_, inertia_in_world);
  }
  throw std::invalid_argument("Error: Unknown frame type. ");
}

void Inertia::checkMass(const double mass) const {
  if (mass <= 0) {
    throw std::invalid_argument("Error: Mass must be positive.");
  }
}

void Inertia::checkInertia(const Eigen::Matrix3d &inertia) const {
  if (!isSymmetric(inertia)) {
    throw std::invalid_argument("Error: Inertia matrix must be symmetric.");
  }
  if (!isPositiveDefinite(inertia)) {
    throw std::invalid_argument("Error: Inertia matrix must be positive semi-definite.");
  }
}

double Inertia::getMass() const { return this->mass_; }
const Eigen::Matrix3d &Inertia::getInertiaTensor() const { return this->inertia_; }
const Frame &Inertia::getFrame() const { return this->frame_; }

geometry_msgs::msg::Inertia Inertia::toRosMessage() const {
  geometry_msgs::msg::Inertia msg;
  msg.m = this->mass_;
  msg.ixx = this->inertia_(0, 0);
  msg.ixy = this->inertia_(0, 1);
  msg.ixz = this->inertia_(0, 2);
  msg.iyy = this->inertia_(1, 1);
  msg.iyz = this->inertia_(1, 2);
  msg.izz = this->inertia_(2, 2);
  return msg;
}
} // namespace spacedyn_ros
