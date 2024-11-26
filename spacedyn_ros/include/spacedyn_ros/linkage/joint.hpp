#ifndef SPACEDYN_ROS_JOINT_HPP_
#define SPACEDYN_ROS_JOINT_HPP_

#include "spacedyn_ros/linkage/joint_state.hpp"
#include "spacedyn_ros/linkage/link_state.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

namespace spacedyn_ros {
class Joint {
public:
  enum class Type { kUndefined = -1, kRevolute = 1, kPrismatic = 2, kFixed = 3 };

  enum ID { kUndefined = -9 };

private:
  int id_;          // Set when connected
  int actuator_id_; // Set when actuator
  std::string name_;
  Transform tf_to_parent_link_com_;   // Only supported local frame expression
  Transform tf_from_parent_link_com_; // Only supported local frame expression

  Type type_;
  Eigen::Vector3d axis_;

  void checkTypeToInput(const Type type) const;

public:
  /***
   * @fn Joint(const Type type)
   * @brief Constructor of Joint class
   * @param type Type of the joint
   * @detail Detailed description
   */
  Joint(const std::string name = "Joint", const Type type = Type::kUndefined,
        const Eigen::Vector3d &axis = Eigen::Vector3d::UnitZ());
  ~Joint() = default;

  // Getters
  std::string getName() const;
  Type getType() const;
  int getId() const;
  int getActuatorId() const;
  bool isActuator() const;

  void connect(const int id, const int actuator_id, const Transform &tf_from_parent_link_com);

  const Transform &getTransformToParentLinkCom() const;
  const Transform &getTransformFromParentLinkCom() const;

  const Eigen::Vector3d &getAxisInLocalFrame() const;
  Eigen::Vector3d getAxisInWorldFrame(const JointState &joint_state) const;
  Eigen::Vector3d getAxisDerivativeInWorldFrame(const JointState &joint_state) const;

  /**
   * @fn transformByActuation(const JointState joint_state)
   * @brief Get the transformation caused by actuation in local frame
   */
  Transform transformByActuation(const JointState &joint_state) const;

  /**
   * @fn twistByActuation(const JointState joint_state)
   * @brief Get the twist caused by actuation in world frame
   */
  Twist twistByActuation(const JointState &joint_state) const;

  /**
   * @fn accelByActuation(const JointState joint_state)
   * @brief Get the acceleration caused by actuation in world frame
   */
  Accel accelByActuation(const JointState &joint_state) const;

  /**
   * @fn computeEffortToAchieveWrench
   * @brief Compute the required effort to achieve the given wrench by the joint.
   */
  double computeEffortToAchieveWrench(const JointState &joint_state) const;

  /**
   * @fn computeVelocityContributionToLinkAccel(const JointState joint_state)
   * @brief Compute the relative twist effect of the joint. This is mainly used to compute the
   * derivative of the jacobian in Kinematics.
   *
   * @param joint_state
   * @return Eigen::Vector6d
   */
  Eigen::Vector6d computeVelocityContributionToLinkAccel(const JointState &joint_state,
                                                         const LinkState &link_state) const;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_JOINT_HPP_
