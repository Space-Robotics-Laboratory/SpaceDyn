#ifndef SPACEDYN_ROS_LINK_HPP_
#define SPACEDYN_ROS_LINK_HPP_

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"
#include "spacedyn_ros/geometry/inertia.hpp"
#include "spacedyn_ros/geometry/transform.hpp"
#include "spacedyn_ros/linkage/joint.hpp"
#include "spacedyn_ros/linkage/link_state.hpp"

namespace spacedyn_ros {
class Link {
private:
  int id_;               // Set when connected
  int parent_id_;        // Set when connected
  int children_number_;  // Update when connected
  bool is_end_effector_; // Update when connected
  std::string name_;
  Inertia inertia_in_local_frame_;

  // Only supported local frame expression
  Transform tf_from_com_to_parent_joint_; // T
  Transform tf_from_parent_joint_to_com_; // T^-1
  Transform tf_com_to_end_tip_;           // Only for end effector

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  enum ID { kUndefined = -9, kBase = 0 };

  /**
   * @fn Link(const std::string name = "Link", const Inertia &inertia_in_local_frame = Inertia(),
   * const Transform &tf_com_to_end_tip = Transform())
   * @brief Construct a new Link object
   * @param name Name of the link
   * @param inertia_in_local_frame Inertia of the link expressed in local frame
   * @param tf_com_to_end_tip Transformation from center of mass to end tip
   * expressed in local frame (only for end effector; otherwise, forced to be an identity matrix)
   */
  Link(const std::string name = "Link", const Inertia &inertia_in_local_frame = Inertia(),
       const Transform &tf_com_to_end_tip = Transform());
  ~Link() = default;

  // Connect
  void connect(const int parent_id, const int assigned_id,
               const Transform &tf_from_parent_joint_to_com);
  void acceptChild();

  // Getters
  std::string getName() const;

  int getId() const;
  int getParentId() const;
  int getParentJointId() const;
  int getChildrenNumber() const;

  double getMass() const;
  const Inertia &getInertiaInLocalFrame() const;
  Inertia computeInertiaInWorldFrame(const Pose &pose) const;

  /**
   * @fn transformParentJoint()
   * @brief Get the transformation matrix from parent joint to the center of
   * mass of the link, expressed in parent joint frame
   */
  const Transform &getTransformToParentJoint() const;
  const Transform &getTransformFromParentJoint() const;
  const Transform &getTransformToEndTip() const;

  bool isBase() const;
  bool isEndEffector() const;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_LINK_HPP_
