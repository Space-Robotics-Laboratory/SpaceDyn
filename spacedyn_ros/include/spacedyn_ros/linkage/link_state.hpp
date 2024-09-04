#ifndef SPACEDYN_ROS_LINK_STATE_HPP_
#define SPACEDYN_ROS_LINK_STATE_HPP_

#include "eigen3/Eigen/Core"
#include "spacedyn_ros/geometry/accel.hpp"
#include "spacedyn_ros/geometry/pose.hpp"
#include "spacedyn_ros/geometry/twist.hpp"
#include "spacedyn_ros/geometry/wrench.hpp"

#include "spacedyn_ros/linkage/joint_state.hpp"

namespace spacedyn_ros {
class LinkState {
private:
  Pose pose_in_world_frame_;                        // position and orientation
  Twist twist_in_world_frame_;                      // velocity of CoM
  Accel accel_in_world_frame_;                      // acceleration of CoM
  Wrench total_wrench_on_link_in_world_frame_;      // summed wrench on link
  Wrench externally_applied_wrench_in_world_frame_; // except gravity

public:
  LinkState();
  ~LinkState() = default;

  const Pose &getPoseInWorldFrame() const;
  const Twist &getTwistInWorldFrame() const;
  const Accel &getAccelInWorldFrame() const;
  const Wrench &getTotalWrenchOnLinkInWorldFrame() const;
  Twist getTwistInLocalFrame() const;
  /**
   * @fn getExternallyAppliedWrenchInWorldFrame()
   * @brief Get the external wrench applied to the link in World frame except for the gravity
   */
  const Wrench &getExternallyAppliedWrenchInWorldFrame() const;

  void setPoseInWorldFrame(const Pose &pose);
  void setTwistInWorldFrame(const Twist &twist);
  void setAccelInWorldFrame(const Accel &accel);
  void setTotalWrenchOnLinkInWorldFrame(const Wrench &wrench);
  void setExternallyAppliedWrenchInWorldFrame(const Wrench &external_wrench);

  void clearPose();
  void clearTwist();
  void clearAccel();
  void clearWrench();
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_LINK_STATE_HPP_
