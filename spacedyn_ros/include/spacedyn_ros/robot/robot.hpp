#ifndef SPACEDYN_ROS_ROBOT_HPP_
#define SPACEDYN_ROS_ROBOT_HPP_

#include "eigen3/Eigen/Core"
#include "geometry_msgs/msg/transform_stamped.hpp"
#include "sensor_msgs/msg/joint_state.hpp"
#include "spacedyn_ros/robot/model.hpp"
#include "spacedyn_ros/robot/state_variable.hpp"

namespace spacedyn_ros {
class Robot { // TODO: Write function description
private:
  Model model_;
  StateVariable state_variable_;

  void checkStateVariable(const StateVariable &state_variable) const;

public:
  Robot(Model const &model = Model());
  Robot(Model const &model, StateVariable const &state_variable);
  Robot(const Robot &robot, const StateVariable &state_variable);
  Robot(const std::string &path_to_urdf);
  ~Robot() = default;

  void operator=(const Robot &robot);

  Eigen::MatrixXd computeGeneralizedJacobianForEndEffector(const int end_effector_id) const;
  Eigen::MatrixXd computeJointToLinkJacobian(const int link_id) const;
  Eigen::MatrixXd computeBaseToLinkJacobian(const int link_id) const;
  Eigen::MatrixXd computeInertiaMatrix() const;
  Eigen::MatrixXd computeGeneralizedInertiaMatrix() const;
  Eigen::VectorXd computeNonlinearVelocityTerm() const;
  Eigen::VectorXd computeGeneralizedNonlinearVelocityTerm() const;
  Eigen::Vector3d computeCenterOfMassInWorldFrame() const;
  Eigen::Vector3d computeVelocityOfCenterOfMassInWorldFrame() const;
  Eigen::Vector3d computeAccelerationOfCenterOfMassInWorldFrame() const;
  Eigen::VectorXd computeMomentumInWorldFrame() const;
  double computeKineticEnergy() const;

  /**
   * @fn computeForward()
   * @brief Compute forward kinematics and dynamics of the robot. The order is, FK(joint pos, vel ->
   * link pose, twist) => FD(link ext wrench -> joint acc) => FK(joint acc -> link acc)
   */
  StateVariable computeForward() const;

  /**
   * @fn step()
   * @brief Step the robot state by one time step. The order is, FK(joint pos, vel -> link pose,
   * twist) => FD(link ext wrench -> joint acc) => FK(joint acc -> link acc) => Integrate(joint pos,
   * vel) => FK(joint pos, vel -> link pose, twist)
   */
  void step();

  const Model &getModel() const;
  double getTotalMass() const;
  const StateVariable &getStateVariable() const;
  int getLinkNumber() const;
  int getJointNumber() const;
  const Link &getLink(const int link_id) const;
  const Link &getEndEffector(const int end_effector_id) const;
  const Joint &getJoint(const int joint_id) const;
  const LinkState &getLinkState(const int link_id) const;
  const JointState &getJointState(const int joint_id) const;
  Eigen::VectorXd getJointPosition() const;
  Eigen::VectorXd getJointVelocity() const;
  Eigen::VectorXd getJointAcceleration() const;
  Eigen::VectorXd getJointEffort() const;
  Eigen::VectorXd getGeneralizedVelocity() const;
  Eigen::VectorXd getGeneralizedAcceleration() const;
  Eigen::VectorXd getGeneralizedForce() const;

  const Pose &getBasePoseInWorldFrame() const;
  const Accel &getBaseAccelInWorldFrame() const;

  const Pose &getLinkPoseInWorldFrame(const int link_id) const;

  // // // Overwrite the state of the robot // // //
  void setStateVariable(const StateVariable &state_variable);
  void overwriteBaseState(const LinkState &base_state);
  void overwriteBasePoseInWorldFrame(const Eigen::Isometry3d &base_pose);
  void overwriteBaseTwistInWorldFrame(const Eigen::VectorXd &base_twist);
  void overwriteBaseAccelInWorldFrame(const Eigen::VectorXd &base_accel);
  void overwriteLinkPoseInWorldFrame(const int link_id, const Eigen::Isometry3d &link_pose);
  void overwriteLinkTwistInWorldFrame(const int link_id, const Eigen::VectorXd &link_twist);
  void overwriteLinkAccelInWorldFrame(const int link_id, const Eigen::VectorXd &link_accel);
  void overwriteJointPosition(const Eigen::VectorXd &joint_position);
  void overwriteJointVelocity(const Eigen::VectorXd &joint_velocity);
  void overwriteJointAcceleration(const Eigen::VectorXd &joint_acceleration);

  void clearBaseState();
  void clearBaseTwist();
  void clearBaseAccel();
  void clearBaseExternallyAppliedWrench();

  void clearLinkState(const int link_id);
  void clearLinkTwist(const int link_id);
  void clearLinkAccel(const int link_id);
  void clearLinkExternallyAppliedWrench(const int link_id);

  void clearAllLinkAccel();
  void clearAllLinkExternallyAppliedWrench();

  void clearJointPosition();
  void clearJointVelocity();
  void clearJointAcceleration();
  void clearJointEffort();

  /**
   * @fn applyExternalWrench()
   * @brief Apply external wrench except gravity to the link
   * @param link_id The id of the link to apply the external wrench
   * @param external_wrench The external wrench to apply (DO NOT include gravity)
   */
  void applyExternalWrench(const int link_id, const Wrench &external_wrench);

  void applyJointEffort(Eigen::VectorXd const &joint_effort);

  // ROS Interface
  geometry_msgs::msg::TransformStamped basePoseToRosTf() const;
  std::vector<geometry_msgs::msg::TransformStamped> jointPoseToRosTf() const;
  void overWriteJointState(const sensor_msgs::msg::JointState &joint_state_msg);
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_ROBOT_HPP_
