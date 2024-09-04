#ifndef SPACEDYN_ROS_DYNAMICS_HPP_
#define SPACEDYN_ROS_DYNAMICS_HPP_

#include "eigen3/Eigen/Core"
#include "spacedyn_ros/geometry/wrench.hpp"
#include "spacedyn_ros/linkage/joint.hpp"
#include "spacedyn_ros/linkage/link.hpp"
#include "spacedyn_ros/robot/robot.hpp"
#include "spacedyn_ros/robot/state_variable.hpp"

namespace spacedyn_ros {
class Dynamics {
private:
  /**
   * @fn inverseLinkWrench
   * @brief Compute the required wrench for the link to achieve the desired acceleration by
   * Newton-Euler equation.
   *
   * @param link
   * @param link_state
   * @return LinkState
   */
  static LinkState inverseLinkWrench(const Link &link, const LinkState &link_state);

  /**
   * @fn inverseJointWrench
   * @brief Compute the required wrench for the joint to satisfy the computed wrench of the link.
   *
   * @param joint
   * @param joint_state
   * @param link
   * @param link_state
   * @param summed_wrench_from_children
   * @param gravity
   * @return JointState
   */
  static JointState inverseJointWrench(const Joint &joint, const JointState &joint_state,
                                       const Link &link, const LinkState &link_state,
                                       const Wrench &summed_wrench_from_children,
                                       const Eigen::Vector3d &gravity);

  /**
   * @fn inverseBaseExtWrench
   * @brief Compute the required external wrench for the base. Used as a substitute for the
   * inverseJointWrench at the base.
   *
   * @param link
   * @param link_state
   * @param summed_wrench_from_children
   * @param gravity
   * @return LinkState
   */
  static LinkState inverseBaseExtWrench(const Link &link, const LinkState &link_state,
                                        const Wrench &summed_wrench_from_children,
                                        const Eigen::Vector3d &gravity);

public:
  Dynamics(/* args */);
  ~Dynamics() = default;

  /**
   * @brief Compute the forward dynamics of the robot. Call Kinematics::computeForward(pose=true,
   * twist=true) first.
   *
   * @param robot
   * @return StateVariable
   */
  static StateVariable computeForward(const Robot &robot);

  /**
   * @fn computeInverse
   * @brief Compute the base wrench and joint effort from the desired link acceleration. Call
   * Kinematics::computeForward(pose=true, twist=true, accel=true) first.
   *
   * @param robot (LinkAccel)
   * @return StateVariable (BaseExternallyAppliedWrench, LinkTotalWrench, JointWrench, JointEffort)
   */
  static StateVariable computeInverse(const Robot &robot);

  // TODO: Check if this function should be here
  static Eigen::Vector3d computeCenterOfMassInWorldFrame(const Robot &robot);
  static Eigen::Vector3d computeVelocityOfCenterOfMassInWorldFrame(const Robot &robot);
  static Eigen::Vector3d computeAccelerationOfCenterOfMassInWorldFrame(const Robot &robot);

  static Eigen::MatrixXd computeInertiaMatrixForBaseMotion(const Robot &robot);
  static Eigen::MatrixXd computeInertiaMatrixForJointMotion(const Robot &robot);
  static Eigen::MatrixXd computeCouplingInertiaMatrix(const Robot &robot);

  static Eigen::MatrixXd computeRobotInertiaMatrix(const Robot &robot);
  static Eigen::MatrixXd computeRobotGeneralizedInertiaMatrix(const Robot &robot);

  static Eigen::VectorXd computeNonlinearVelocityTerm(const Robot &robot);
  static Eigen::VectorXd computeGeneralizedNonlinearVelocityTerm(const Robot &robot);

  static Eigen::VectorXd computeRobotMomentumInWorldFrame(const Robot &robot);

  static double computeRobotKineticEnergy(const Robot &robot);
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_DYNAMICS_HPP_
