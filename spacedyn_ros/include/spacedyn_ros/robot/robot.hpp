#ifndef SPACEDYN_ROS_ROBOT_HPP_
#define SPACEDYN_ROS_ROBOT_HPP_

#include "spacedyn_ros/robot/model.hpp"
#include "spacedyn_ros/robot/state_variable.hpp"
#include <eigen3/Eigen/Core>
#include <geometry_msgs/msg/transform_stamped.hpp>
#include <sensor_msgs/msg/joint_state.hpp>

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

  /**
   * @fn computeGeneralizedJacobianForLink()
   * @brief Return GJ (6, n) = J_mi - J_bi * H_b^-1 * H_bmi
   * @param link_id The id of the link to compute the generalized jacobian
   * @return Eigen::MatrixXd GJ (6, n)
   */
  Eigen::MatrixXd computeGeneralizedJacobianForLink(const int link_id) const;
  /**
   * @fn computeGeneralizedJacobianForEndEffector()
   * @brief Return GJ (6, n) = J_mi - J_bi * H_b^-1 * H_bmi
   * @param end_effector_id The id of the End-Effector to compute the generalized jacobian
   * @return Eigen::MatrixXd GJ (6, n)
   */
  Eigen::MatrixXd computeGeneralizedJacobianForEndEffector(const int end_effector_id) const;
  /**
   * @fn computeGeneralizedJacobianForEndTip()
   * @brief Return GJ (6, n) = J_mi - J_bi * H_b^-1 * H_bmi
   * @param end_effector_id The id of the End-Effector to compute the generalized jacobian
   * @return Eigen::MatrixXd GJ (6, n)
   */
  Eigen::MatrixXd computeGeneralizedJacobianForEndTip(const int end_effector_id) const;
  /**
   * @fn computeJointToLinkJacobian()
   * @brief Compute the Jacobian matrix of the end effector effected by joints motion with respect
   * to the base frame. Call this function after calling computeForward() with compute_pose = true.
   * @param link_id The id of the link to compute jacobian
   * @return Eigen::MatrixXd J_mi (6, n)
   */
  Eigen::MatrixXd computeJointToLinkJacobian(const int link_id) const;
  /**
   * @fn computeJointToEndTipJacobian()
   * @brief Compute the Jacobian matrix of the tip point of the end effector effected by joints
   * motion with respect to the base frame. End tip is set when the end-effector link was generated.
   * Call this function after calling computeForward() with compute_pose = true.
   * @param end_effector_id The id of the End-Effector to compute jacobian
   * @return Eigen::MatrixXd J_mi (6, n)
   */
  Eigen::MatrixXd computeJointToEndTipJacobian(const int end_effector_id) const;
  /**
   * @fn computeBaseToLinkJacobian()
   * @brief Compute the Jacobian matrix of the end effector effected by the base motion with respect
   * to the base frame. Call this function after calling computeForward() with compute_pose = true.
   * @param link_id The id of the link to compute jacobian
   * @return Eigen::MatrixXd J_bi (6, n)
   */
  Eigen::MatrixXd computeBaseToLinkJacobian(const int link_id) const;
  /**
   * @fn computeBaseToEndTipJacobian()
   * @brief Compute the Jacobian matrix of the tip point of the end effector effected by the base
   * motion with respect to the base frame. End tip is set when the end-effector link was generated.
   * Call this function after calling computeForward() with compute_pose = true.
   * @param end_effector_id The id of the End-Effector to compute jacobian
   * @return Eigen::MatrixXd J_bi (6, n)
   */
  Eigen::MatrixXd computeBaseToEndTipJacobian(const int end_effector_id) const;
  /**
   * @fn computeInertiaMatrix()
   * @brief Compute the robot inertia matrix.
   * @return Eigen::MatrixXd H (6+n, 6+n) = [H_b, H_bm; H_bm^T, H_m]
   */
  Eigen::MatrixXd computeInertiaMatrix() const;
  /**
   * @fn computeGeneralizedInertiaMatrix()
   * @brief Compute the generalized inertia matrix.
   * @return Eigen::MatrixXd GH (n, n) = H_m - H_bm^T * H_b^-1 * H_bm
   */
  Eigen::MatrixXd computeGeneralizedInertiaMatrix() const;
  Eigen::MatrixXd computeInertiaMatrixForBaseMotion() const;
  Eigen::MatrixXd computeCouplingInertiaMatrix() const;
  Eigen::MatrixXd computeInertiaMatrixForJointMotion() const;
  /**
   * @fn computeNonlinearVelocityTerm()
   * @brief Compute the nonlinear velocity term.
   * @return Eigen::VectorXd C (6+n, 1) = [C_b; C_m];
   */
  Eigen::VectorXd computeNonlinearVelocityTerm() const;
  /**
   * @fn computeGeneralizedNonlinearVelocityTerm()
   * @brief Compute the generalized nonlinear velocity term.
   * @return Eigen::VectorXd GC (n, n) = C_m - H_bm^T * H_b^-1 * C_b
   */
  Eigen::VectorXd computeGeneralizedNonlinearVelocityTerm() const;
  /**
   * @fn computeGravityTerm()
   * @brief Compute the gravity term.
   * @return Eigen::VectorXd G (6+n, 1)
   */
  Eigen::VectorXd computeGravityTerm() const;
  /**
   * @fn computeCenterOfMassInWorldFrame()
   * @brief Compute the center of mass of the robot in the world frame.
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d computeCenterOfMassInWorldFrame() const;
  /**
   * @fn computeVelocityOfCenterOfMassInWorldFrame()
   * @brief Compute the velocity of the center of mass of the robot in the world frame.
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d computeVelocityOfCenterOfMassInWorldFrame() const;
  /**
   * @fn computeAccelerationOfCenterOfMassInWorldFrame()
   * @brief Compute the acceleration of the center of mass of the robot in the world frame.
   * @return Eigen::Vector3d
   */
  Eigen::Vector3d computeAccelerationOfCenterOfMassInWorldFrame() const;
  /**
   * @fn computeMomentumInWorldFrame()
   * @brief Compute the momentum of the robot in the world frame.
   * @return Eigen::Vector6d
   */
  Eigen::Vector6d computeMomentumInWorldFrame() const;
  Eigen::Vector6d computeMomentumAroundBaseInWorldFrame() const;
  /**
   * @fn computeKineticEnergy()
   * @brief Compute the kinetic energy of the robot.
   * @return double
   */
  double computeKineticEnergy() const;

  /**
   * @fn computeForward()
   * @brief Compute forward kinematics and dynamics of the robot. The order is, FK(joint pos, vel ->
   * link pose, twist) => FD(link ext wrench -> joint acc) => FK(joint acc -> link acc)
   */
  StateVariable computeForward() const;

  /**
   * @fn computeInverseDynamics()
   * @brief Compute inverse dynamics of the robot.
   * @return VectorXd(6+n, 1) = [BaseExternallyAppliedWrench; JointEffort]
   */
  Eigen::VectorXd computeInverseDynamics(const Eigen::Vector6d &desired_base_acceleration,
                                         const Eigen::VectorXd &desired_joint_acceleration) const;

  /**
   * @fn computeInverseDynamicsInJointSpace()
   * @brief Compute inverse dynamics of the robot in joint space.
   * @return VectorXd(n, 1) = JointEffort
   */
  Eigen::VectorXd
  computeInverseDynamicsInJointSpace(const Eigen::VectorXd &desired_joint_acceleration) const;

  const Model &getModel() const;
  double getTotalMass() const;
  const StateVariable &getStateVariable() const;

  /**
   * @fn getDof()
   * @brief Return the degree of freedom of the robot: 6 (base) + n (actuator)
   * @return int
   */
  int getDof() const;
  int getLinkNumber() const;
  int getJointNumber() const;
  int getActuatorNumber() const;
  int getEndEffectorNumber() const;
  const Link &getLink(const int link_id) const;
  const Link &getEndEffector(const int end_effector_id) const;
  const Joint &getJoint(const int joint_id) const;
  const LinkState &getLinkState(const int link_id) const;
  const LinkState &getEndEffectorState(const int end_effector_id) const;
  const JointState &getJointState(const int joint_id) const;
  Eigen::VectorXd getJointPosition() const;
  Eigen::VectorXd getJointVelocity() const;
  Eigen::VectorXd getJointAcceleration() const;
  Eigen::VectorXd getJointEffort() const;

  /**
   * @fn getGeneralizedVelocity()
   * @brief Return GV (6+n, 1) = [base_twist; joint_velocity]
   */
  Eigen::VectorXd getGeneralizedVelocity() const;

  /**
   * @fn getGeneralizedAcceleration()
   * @brief Return GA (6+n, 1) = [base_accel; joint_acceleration]
   */
  Eigen::VectorXd getGeneralizedAcceleration() const;

  /**
   * @fn getGeneralizedForce()
   * @brief Return GF (6+n, 1) = [base_wrench; joint_effort] + sum ([J_bi^T, J_mi^T] * link_wrench)
   */
  Eigen::VectorXd getGeneralizedForce() const;

  const Eigen::Isometry3d &getBasePose() const; // TODO: Remove this function
  const Pose &getBasePoseInWorldFrame() const;
  const Pose &getLinkPoseInWorldFrame(const int link_id) const;
  const Pose &getEndEffectorPoseInWorldFrame(const int end_effector_id) const;
  const Pose getEndTipPoseInWorldFrame(const int end_effector_id) const;

  const Twist &getBaseTwistInWorldFrame() const;
  Twist getBaseTwistInLocalFrame() const;
  const Twist &getLinkTwistInWorldFrame(const int link_id) const;
  Twist getLinkTwistInLocalFrame(const int link_id) const;
  const Twist &getEndEffectorTwistInWorldFrame(const int end_effector_id) const;
  Twist getEndEffectorTwistInLocalFrame(const int end_effector_id) const;
  Twist getEndTipTwistInWorldFrame(const int end_effector_id) const;
  // TODO: Add getEndTipTwistInLocalFrame()

  const Accel &getBaseAccelInWorldFrame() const;
  Accel getBaseAccelInLocalFrame() const;
  const Accel &getLinkAccelInWorldFrame(const int link_id) const;
  Accel getLinkAccelInLocalFrame(const int link_id) const;
  const Accel &getEndEffectorAccelInWorldFrame(const int end_effector_id) const;
  Accel getEndEffectorAccelInLocalFrame(const int end_effector_id) const;
  Accel getEndTipAccelInWorldFrame(const int end_effector_id) const;
  // TODO: Add getEndTipAccelInLocalFrame()

  // // // Overwrite the state of the robot // // //
  /**
   * @fn step()
   * @brief Step the robot state by one time step. The order is, FK(joint pos, vel -> link pose,
   * twist) => FD(link ext wrench -> joint acc) => FK(joint acc -> link acc) => Integrate(joint pos,
   * vel) => FK(joint pos, vel -> link pose, twist)
   */
  void step();
  void updateKinematics(bool pose = true, bool twist = true, bool accel = true);
  void setStateVariable(const StateVariable &state_variable);
  void overwriteBaseState(const LinkState &base_state);
  void overwriteBasePoseInWorldFrame(const Eigen::Isometry3d &base_pose);
  void overwriteBasePoseInWorldFrame(const Eigen::Vector3d &base_position,
                                     const Eigen::Matrix3d &base_attitude);
  void overwriteBaseTwistInWorldFrame(const Eigen::Vector6d &base_twist);
  void overwriteBaseAccelInWorldFrame(const Eigen::Vector6d &base_accel);
  void overwriteLinkPoseInWorldFrame(const int link_id, const Eigen::Isometry3d &link_pose);
  void overwriteLinkTwistInWorldFrame(const int link_id, const Eigen::Vector6d &link_twist);
  void overwriteLinkAccelInWorldFrame(const int link_id, const Eigen::Vector6d &link_accel);
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

  /**
   * @fn applyJointEffort()
   * @brief Apply joint effort exerted by the joint actuators
   * @param joint_effort The joint effort to apply
   */
  void applyJointEffort(Eigen::VectorXd const &joint_effort);

  // ROS Interface
  geometry_msgs::msg::TransformStamped basePoseToRosTf() const;
  std::vector<geometry_msgs::msg::TransformStamped> jointPoseToRosTf() const;
  void overWriteJointState(const sensor_msgs::msg::JointState &joint_state_msg);
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_ROBOT_HPP_
