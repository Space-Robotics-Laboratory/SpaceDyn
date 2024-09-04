#include "spacedyn_ros/linkage/joint_state.hpp"
#include "eigen3/Eigen/Core"
#include "iostream"
#include "spacedyn_ros/geometry/accel.hpp"
#include "spacedyn_ros/geometry/pose.hpp"
#include "spacedyn_ros/geometry/twist.hpp"

namespace spacedyn_ros {

JointState::JointState() {
  // Initialize the joint state with the given joint number
  this->position_ = 0;
  this->velocity_ = 0;
  this->acceleration_ = 0;
  this->effort_ = 0;

  this->pose_in_world_frame_ = Pose();
  this->twist_in_world_frame_ = Twist();
  this->accel_in_world_frame_ = Accel();
  this->wrench_to_child_in_world_frame_ = Wrench();
}

const Pose &JointState::getPoseInWorldFrame() const { return this->pose_in_world_frame_; }
const Twist &JointState::getTwistInWorldFrame() const { return this->twist_in_world_frame_; }
const Accel &JointState::getAccelInWorldFrame() const { return this->accel_in_world_frame_; }
const Wrench &JointState::getWrenchToChildInWorldFrame() const {
  return this->wrench_to_child_in_world_frame_;
}

double JointState::getPosition() const { return this->position_; }
double JointState::getVelocity() const { return this->velocity_; }
double JointState::getAcceleration() const { return this->acceleration_; }
double JointState::getEffort() const { return this->effort_; }

Eigen::Vector3d JointState::getAxisInWorldFrame() const {
  return this->pose_in_world_frame_.getOriginPose().linear().col(2);
}

Eigen::Vector3d JointState::getAxisDerivativeInWorldFrame() const {
  return this->twist_in_world_frame_.getOriginAngularVelocity().cross(this->getAxisInWorldFrame());
}

void JointState::setPoseInWorldFrame(const Pose &pose) { this->pose_in_world_frame_ = pose; }
void JointState::setTwistInWorldFrame(const Twist &twist) {
  if (twist.getFrame() != Frame::kWorld) {
    throw std::invalid_argument("Error: twist should be in world frame. ");
  }
  this->twist_in_world_frame_ = twist;
}
void JointState::setAccelInWorldFrame(const Accel &accel) {
  if (accel.getFrame() != Frame::kWorld) {
    throw std::invalid_argument("Error: accel should be in world frame. ");
  }
  this->accel_in_world_frame_ = accel;
}

void JointState::setWrenchToChildInWorldFrame(const Wrench &wrench) {
  if (wrench.getFrame() != Frame::kWorld) {
    throw std::invalid_argument("Error: wrench should be in world frame. ");
  }
  this->wrench_to_child_in_world_frame_ = wrench;
}

void JointState::setPosition(double position) { this->position_ = position; }
void JointState::setVelocity(double velocity) { this->velocity_ = velocity; }
void JointState::setAcceleration(double acceleration) { this->acceleration_ = acceleration; }
void JointState::setEffort(double effort) { this->effort_ = effort; }

void JointState::clearPose() { this->pose_in_world_frame_ = Pose(); }
void JointState::clearTwist() { this->twist_in_world_frame_ = Twist(); }
void JointState::clearAccel() { this->accel_in_world_frame_ = Accel(); }

} // namespace spacedyn_ros
