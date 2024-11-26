#ifndef SPACEDYN_ROS_KINEMATICS_HPP_
#define SPACEDYN_ROS_KINEMATICS_HPP_

#include "spacedyn_ros/linkage/link_state.hpp"
#include "spacedyn_ros/robot/robot.hpp"
#include "spacedyn_ros/robot/state_variable.hpp"
#include <eigen3/Eigen/Core>

namespace spacedyn_ros {
class Kinematics {
private:
  /**
   * @fn JointState forwardJointPose
   * @brief This function computes the forward kinematics about Pose of the
   * given joint.
   * @param (joint) Joint model whose pose is computed
   * @param (parent_link_state) Link state of the parent link
   * @param (joint_state) Joint state whose pose is computed
   * @return Returns joint state with the updated pose. Twist and accel are not
   * updated.
   */
  static JointState forwardJointPose(const Joint &joint, const LinkState &parent_link_state,
                                     const JointState &joint_state);
  /**
   * @fn JointState forwardJointTwist
   * @brief This function computes the forward kinematics about Twist of the
   * given joint.
   * @param (joint) Joint model whose twist is computed
   * @param (parent_link_state) Link state of the parent link
   * @param (joint_state) Joint state whose twist is computed
   * @return Returns joint state with the updated twist. Pose and accel are not
   * updated.
   */
  static JointState forwardJointTwist(const Joint &joint, const LinkState &parent_link_state,
                                      const JointState &joint_state);
  /**
   * @fn JointState forwardJointAccel
   * @brief This function computes the forward kinematics about Accel of the
   * given joint.
   * @param (joint) Joint model whose accel is computed
   * @param (parent_link_state) Link state of the parent link
   * @param (joint_state) Joint state whose accel is computed
   * @return Returns joint state with the updated accel. Pose and twist are not
   * updated.
   */
  static JointState forwardJointAccel(const Joint &joint, const LinkState &parent_link_state,
                                      const JointState &joint_state);

  /**
   * @fn LinkState forwardLinkPose(Link link, JointState parent_joint_state,
   * LinkState link_state)
   * @brief This function computes the forward kinematics about Pose of the
   * given link. Call forwardJointPose() first.
   * @param (link) Link model whose pose is computed
   * @param (parent_joint_state) Joint state of the parent joint
   * @param (link_state) Link state whose pose is computed
   * @return Returns link state with the updated pose. Twist and accel are not
   * updated.
   */
  static LinkState forwardLinkPose(const Link &link, const JointState &parent_joint_state,
                                   const LinkState &link_state);

  /**
   * @fn LinkState forwardLinkTwist
   * @brief This function computes the forward kinematics about Twist of the
   * given link. Call forwardJointTwist() first.
   * @param (link) Link model whose twist is computed
   * @param (parent_joint_state) Joint state of the parent joint
   * @param (link_state) Link state whose twist is computed
   * @return Returns link state with the updated twist. Pose and accel are not
   * updated.
   */
  static LinkState forwardLinkTwist(const Link &link, const JointState &parent_joint_state,
                                    const LinkState &link_state);
  /**
   * @fn LinkState forwardLinkAccel
   * @brief This function computes the forward kinematics about Accel of the
   * given link. Call forwardJointAccel() first.
   * @param (link) Link model whose accel is computed
   * @param (parent_joint_state) Joint state of the parent joint
   * @param (link_state) Link state whose accel is computed
   * @return Returns link state with the updated accel. Pose and twist are not
   * updated.
   */
  static LinkState forwardLinkAccel(const Link &link, const JointState &parent_joint_state,
                                    const LinkState &link_state);

public:
  Kinematics(/* args */);
  ~Kinematics() = default;

  /**
   * @fn computeForward
   * @brief Compute the forward kinematics of the robot. This function computes the pose, twist, and
   * accel of each link and joint recursively in the robot. The result is stored in the
   * StateVariable object. The computation is done in the order of pose -> twist -> accel.
   *
   * @param robot
   * @param compute_pose
   * @param compute_twist
   * @param compute_accel
   * @return StateVariable
   */
  static StateVariable computeForward(const Robot &robot, const bool compute_pose,
                                      const bool compute_twist, const bool compute_accel);
  static StateVariable computeInverse(const Robot &robot);

  static Eigen::MatrixXd computeGeneralizedJacobianForLink(const Robot &robot, const int link_id);
  static Eigen::MatrixXd computeGeneralizedJacobianForEndEffector(const Robot &robot,
                                                                  const int end_effector_id);
  static Eigen::MatrixXd computeGeneralizedJacobianForEndTip(const Robot &robot,
                                                             const int end_effector_id);

  // TODO: Check if the name is suitable
  /**
   * @brief Compute the Jacobian matrix of the end effector effected by joints motion with respect
   * to the base frame. Call this function after calling computeForward() with compute_pose = true.
   *
   * @param robot
   * @param link_id
   * @return Eigen::MatrixXd
   */
  static Eigen::MatrixXd computeJointToLinkJacobian(const Robot &robot, const int link_id);

  /**
   * @fn computeJointToEndTipJacobian
   * @brief Compute the Jacobian matrix of the tip point of the end effector effected by joints
   * motion with respect to the base frame. End tip is set when the end-effector link was generated.
   * Call this function after calling computeForward() with compute_pose = true.
   */
  static Eigen::MatrixXd computeJointToEndTipJacobian(const Robot &robot,
                                                      const int end_effector_id);

  /**
   * @fn computeJointToLinkJacobianDerivative
   * @brief Compute the derivative of the Jacobian matrix of the end effector effected by joints
   * motion with respect to the base frame. Call this function after calling computeForward() with
   * compute_pose = true, compute_twist = true.
   *
   * @param robot
   * @param link_id
   * @return Eigen::MatrixXd
   */
  static Eigen::MatrixXd computeJointToLinkJacobianDerivative(const Robot &robot,
                                                              const int link_id);

  /**
   * @brief Compute the Jacobian matrix of the end effector effected by the base motion with respect
   * to the base frame. Call this function after calling computeForward() with compute_pose = true.
   *
   * @param robot
   * @param link_id
   * @return Eigen::MatrixXd
   */
  static Eigen::MatrixXd computeBaseToLinkJacobian(const Robot &robot, const int link_id);

  /**
   * @fn computeBaseToEndTipJacobian
   * @brief Compute the Jacobian matrix of the tip point of the end effector effected by the base
   * motion with respect to the base frame. End tip is set when the end-effector link was generated.
   * Call this function after calling computeForward() with compute_pose = true.
   */
  static Eigen::MatrixXd computeBaseToEndTipJacobian(const Robot &robot, const int end_effector_id);

  /**
   * @fn computeBaseToLinkJacobianDerivative
   * @brief Compute the derivative of the Jacobian matrix of the end effector effected by the base
   * motion with respect to the base frame. Call this function after calling computeForward() with
   * compute_pose = true, compute_twist = true.
   *
   * @param robot
   * @param link_id
   * @return Eigen::MatrixXd
   */
  static Eigen::MatrixXd computeBaseToLinkJacobianDerivative(const Robot &robot, const int link_id);
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_KINEMATICS_HPP_
