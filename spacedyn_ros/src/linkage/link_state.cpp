#include "spacedyn_ros/linkage/link_state.hpp"
#include "eigen3/Eigen/Core"

#include "iostream"

namespace spacedyn_ros {

LinkState::LinkState() {
  this->pose_in_world_frame_ = Pose();
  this->twist_in_world_frame_ = Twist();
  this->accel_in_world_frame_ = Accel();
  this->total_wrench_on_link_in_world_frame_ = Wrench();
}

const Pose &LinkState::getPoseInWorldFrame() const { return this->pose_in_world_frame_; }
const Twist &LinkState::getTwistInWorldFrame() const { return this->twist_in_world_frame_; }
Twist LinkState::getTwistInLocalFrame() const {
  auto twist = getPoseInWorldFrame().computeTwistInLocalFrame(getTwistInWorldFrame());
  return twist;
}
const Accel &LinkState::getAccelInWorldFrame() const { return this->accel_in_world_frame_; }
const Wrench &LinkState::getTotalWrenchOnLinkInWorldFrame() const {
  return this->total_wrench_on_link_in_world_frame_;
}
const Wrench &LinkState::getExternallyAppliedWrenchInWorldFrame() const {
  return this->externally_applied_wrench_in_world_frame_;
}

void LinkState::setPoseInWorldFrame(const Pose &pose) { this->pose_in_world_frame_ = pose; }
void LinkState::setTwistInWorldFrame(const Twist &twist) {
  if (twist.getFrame() != Frame::kWorld) {
    throw std::invalid_argument("Error: twist should be in world frame. ");
  }
  this->twist_in_world_frame_ = twist;
}
void LinkState::setAccelInWorldFrame(const Accel &accel) {
  if (accel.getFrame() != Frame::kWorld) {
    throw std::invalid_argument("Error: accel should be in world frame. ");
  }
  this->accel_in_world_frame_ = accel;
}
void LinkState::setTotalWrenchOnLinkInWorldFrame(const Wrench &wrench) {
  if (wrench.getFrame() != Frame::kWorld) {
    throw std::invalid_argument("Error: wrench should be in world frame. ");
  }
  this->total_wrench_on_link_in_world_frame_ = wrench;
}

void LinkState::setExternallyAppliedWrenchInWorldFrame(const Wrench &external_wrench) {
  if (external_wrench.getFrame() == Frame::kLocal) {
    this->externally_applied_wrench_in_world_frame_ =
        pose_in_world_frame_.computeWrenchInWorldFrame(external_wrench);
  } else if (external_wrench.getFrame() == Frame::kWorld) {
    this->externally_applied_wrench_in_world_frame_ = external_wrench;
  }
}

void LinkState::clearPose() { this->pose_in_world_frame_ = Pose(); }
void LinkState::clearTwist() { this->twist_in_world_frame_ = Twist(); }
void LinkState::clearAccel() { this->accel_in_world_frame_ = Accel(); }
void LinkState::clearWrench() { this->total_wrench_on_link_in_world_frame_ = Wrench(); }

} // namespace spacedyn_ros
