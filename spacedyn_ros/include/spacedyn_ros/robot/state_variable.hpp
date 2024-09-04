#ifndef SPACEDYN_ROS_LINK_STATE_VARIABLE_HPP_
#define SPACEDYN_ROS_LINK_STATE_VARIABLE_HPP_

#include "eigen3/Eigen/Core"
#include "spacedyn_ros/linkage/joint_state.hpp"
#include "spacedyn_ros/linkage/link_state.hpp"
#include "spacedyn_ros/robot/model.hpp"

// TODO: Add variable to handle state update condition

namespace spacedyn_ros {
class StateVariable {
private:
  std::vector<LinkState> links_state_;
  std::vector<JointState> joints_state_;
  int joint_number_;
  int link_number_;

  void checkModelToSet(const Model &model) const;
  void checkLinkIdToCall(const int link_id) const;
  void checkJointIdToCall(const int joint_id) const;
  void checkJointVectorToSet(const Eigen::VectorXd &joint_vector) const;
  void checkGeneralizedVectorToSet(const Eigen::VectorXd &generalized_vector) const;

public:
  StateVariable(Model const &model);
  StateVariable(const int link_number);
  ~StateVariable() = default;

  StateVariable copyVacant() const;

  int getLinkNumber() const;
  int getJointNumber() const;

  const LinkState &getLinkState(const int link_id) const;
  const LinkState &getLinkState(const Link &link) const;
  const JointState &getJointState(const int joint_id) const;
  const JointState &getJointState(const Joint &joint) const;

  // Handle all joints' position, velocity, acceleration, and effort as
  // Eigen::VectorXd at the same time
  Eigen::VectorXd getJointPosition() const;
  Eigen::VectorXd getJointVelocity() const;
  Eigen::VectorXd getJointAcceleration() const;
  Eigen::VectorXd getJointEffort() const;

  Eigen::VectorXd getGeneralizedVelocity() const;
  Eigen::VectorXd getGeneralizedAcceleration() const;

  void setLinkPoseInWorldFrame(const int link_id, const Eigen::Isometry3d &pose);
  void setLinkTwistInWorldFrame(const int link_id, const Eigen::VectorXd &twist);
  void setLinkAccelInWorldFrame(const int link_id, const Eigen::VectorXd &accel);
  void setLinkExternallyAppliedWrenchInWorldFrame(const int link_id, const Eigen::VectorXd &wrench);

  void setJointPosition(const Eigen::VectorXd &joint_position);
  void setJointVelocity(const Eigen::VectorXd &joint_velocity);
  void setJointAcceleration(const Eigen::VectorXd &joint_acceleration);
  void setJointEffort(const Eigen::VectorXd &joint_effort);

  void setLinkState(const int link_id, const LinkState &link_state);
  void setJointState(const int joint_id, const JointState &joint_state);

  void setGeneralizedCoordinates(const Eigen::Isometry3d &base_pose,
                                 const Eigen::VectorXd &joint_position);
  void setGeneralizedVelocity(const Eigen::VectorXd &generalized_velocity);
  void setGeneralizedAcceleration(const Eigen::VectorXd &generalized_acceleration);
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_LINK_STATE_VARIABLE_HPP_
