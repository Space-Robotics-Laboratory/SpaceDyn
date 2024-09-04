#ifndef SPACEDYN_ROS_JOINT_STATE_HPP_
#define SPACEDYN_ROS_JOINT_STATE_HPP_

#include "eigen3/Eigen/Core"
#include "spacedyn_ros/geometry/accel.hpp"
#include "spacedyn_ros/geometry/pose.hpp"
#include "spacedyn_ros/geometry/twist.hpp"
#include "spacedyn_ros/geometry/wrench.hpp"

namespace spacedyn_ros {
class JointState {
private:
  // Joint state as rigid local
  // Frame is fixed to child link
  Pose pose_in_world_frame_;
  Twist twist_in_world_frame_;
  Accel accel_in_world_frame_;
  Wrench wrench_to_child_in_world_frame_;

  // Joint state as actuator
  double position_;
  double velocity_;
  double acceleration_;
  double effort_;

public:
  JointState();
  ~JointState() = default;

  /**
   * @fn getPose()
   * @brief Get the pose of the joint after actuation in World frame
   */
  const Pose &getPoseInWorldFrame() const;
  const Twist &getTwistInWorldFrame() const;
  const Accel &getAccelInWorldFrame() const;
  const Wrench &getWrenchToChildInWorldFrame() const;

  Eigen::Vector3d getAxisInWorldFrame() const;
  Eigen::Vector3d getAxisDerivativeInWorldFrame() const;

  double getPosition() const;
  double getVelocity() const;
  double getAcceleration() const;
  double getEffort() const;

  void setPoseInWorldFrame(const Pose &pose);
  void setTwistInWorldFrame(const Twist &twist);
  void setAccelInWorldFrame(const Accel &accel);
  void setWrenchToChildInWorldFrame(const Wrench &wrench);

  void setPosition(double position);
  void setVelocity(double velocity);
  void setAcceleration(double acceleration);
  void setEffort(double effort);

  void clearPose();
  void clearTwist();
  void clearAccel();
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_JOINT_STATE_HPP_
