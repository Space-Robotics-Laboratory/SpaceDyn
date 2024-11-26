#ifndef SPACEDYN_ROS_DYNAMICS_HPP_
#define SPACEDYN_ROS_DYNAMICS_HPP_

#include "spacedyn_ros/geometry/wrench.hpp"
#include "spacedyn_ros/linkage/joint.hpp"
#include "spacedyn_ros/linkage/link.hpp"
#include "spacedyn_ros/robot/robot.hpp"
#include "spacedyn_ros/robot/state_variable.hpp"
#include <eigen3/Eigen/Core>

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
   * @fn computeForward
   * @brief Compute the link acceleration from the base wrench and joint effort.
   * Compute the q_ddot = H^-1 * (Ta - C(q, q_dot) - G(q))
   *
   * @param robot (LinkExternallyAppliedWrench, JointEffort)
   * @return StateVariable (BaseAccel, joint_accel)
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

  /**
   * @fn computeInverseInJointSpace
   * @brief Compute the joint effort from the desired joint acceleration. Call
   * Kinematics::computeForward(pose=true, twist=true, accel=true) first.
   */
  static Eigen::VectorXd computeInverseInJointSpace(const Robot &robot);

  // TODO: Check if this function should be here
  /**
   * @fn computeCenterOfMassInWorldFrame
   * @brief Compute the center of mass of the robot in the world frame.
   */
  static Eigen::Vector3d computeCenterOfMassInWorldFrame(const Robot &robot);
  static Eigen::Vector3d computeVelocityOfCenterOfMassInWorldFrame(const Robot &robot);
  static Eigen::Vector3d computeAccelerationOfCenterOfMassInWorldFrame(const Robot &robot);

  /**
   * @fn computeInertiaMatrixForBaseMotion
   * @brief Compute the base part of the robot inertia matrix for the base motion. This is used
   * to compute the effect of the base motion as a part of the robot system.
   *
   * @param robot
   * @return Hb (6, 6)
   */
  static Eigen::MatrixXd computeInertiaMatrixForBaseMotion(const Robot &robot);

  /**
   * @fn computeInertiaMatrixForJointMotion
   * @brief Compute the joint part of the robot inertia matrix for the joint motion. This is used to
   * compute the effect of the joint motion on the joint motion as a part of the robot system.
   *
   * @param robot
   * @return Hm (n, n)
   */
  static Eigen::MatrixXd computeInertiaMatrixForJointMotion(const Robot &robot);

  /**
   * @fn computeCouplingInertiaMatrix
   * @brief Compute the coupling inertia matrix. This is used to compute the interaction between the
   * base and the joint motion.
   *
   * @param robot
   * @return Hbm (6, n)
   */
  static Eigen::MatrixXd computeCouplingInertiaMatrix(const Robot &robot);

  /**
   * @fn computeRobotInertiaMatrix
   * @brief Compute the robot inertia matrix.
   *
   * @param robot
   * @return H (6+n, 6+n) = [H_b, H_bm; H_bm^T, H_m].
   */
  static Eigen::MatrixXd computeRobotInertiaMatrix(const Robot &robot);
  static Eigen::MatrixXd computeRobotGeneralizedInertiaMatrix(const Robot &robot);

  /**
   * @fn computeNonlinearVelocityTerm
   * @brief Compute the nonlinear velocity term of the robot dynamics using recursive Newton-Euler
   * method.
   *
   * @param robot
   * @return C(q, q_dot) (6+n, 1)
   */
  static Eigen::VectorXd computeNonlinearVelocityTerm(const Robot &robot);
  static Eigen::VectorXd computeGeneralizedNonlinearVelocityTerm(const Robot &robot);

  /**
   * @fn computeGravityTerm
   * @brief Compute the gravity term of the robot dynamics.
   *
   * @param robot
   * @return G(q) (6+n, 1)
   */
  static Eigen::VectorXd computeGravityTerm(const Robot &robot);
  static Eigen::VectorXd computeGeneralizedGravityTerm(const Robot &robot);

  static Eigen::Vector6d computeRobotMomentumInWorldFrame(const Robot &robot);
  static Eigen::Vector6d computeRobotMomentumAroundBaseInWorldFrame(const Robot &robot);

  /**
   * @fn computeRobotKineticEnergy
   * @brief Compute the kinetic energy of the robot.
   *
   * @param robot
   * @return double
   */
  static double computeRobotKineticEnergy(const Robot &robot);
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_DYNAMICS_HPP_
