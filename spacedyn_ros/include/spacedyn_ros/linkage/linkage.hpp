#ifndef SPACEDYN_ROS_LINKAGE_HPP_
#define SPACEDYN_ROS_LINKAGE_HPP_

#include "eigen3/Eigen/Core"
#include "spacedyn_ros/linkage/link.hpp"
#include <vector>

namespace spacedyn_ros {
class Linkage {
private:
  // Joint and Link object
  std::vector<Joint> joints_;
  std::vector<Link> links_;

  // Physical properties
  double total_mass_;

  // ID for component
  std::vector<int> end_effectors_;
  std::vector<std::vector<int>> link_id_jag_to_end_effector_;

  // ID checker
  void checkLinkIdToCall(const int link_id) const;
  void checkJointIdToCall(const int joint_id) const;
  void checkEndEffectorIdToCall(const int end_effector_id) const;

  // Add link and joint to the linkage

  void addLink(const int parent_link_id, const Link &child_input,
               const Transform &tf_from_parent_joint_to_com);

  void addJoint(const int parent_link_id, const Joint &joint_input,
                const Transform &tf_from_parent_link_com_to_joint);
  void replaceLink(const int id, const Link &link);

  std::vector<int> computeLinkIdChainToEndEffector(const int end_effector_id) const;

  std::vector<int> replaceEndEffector(const Link parent, const Link child) const;

public:
  Linkage(/* args */);
  ~Linkage() = default;

  void addBase(const Link &link);

  /**
   * @fn void addLink(const int parent_link_id, const Link &child_input, const Transform
   * &tf_from_parent_joint_to_com)
   * @brief Add link to linkage.
   */
  void addJointWithLink(const int parent_link_id, const Joint &child_joint,
                        const Transform &tf_from_parent_link_com, const Link &child_link,
                        const Transform &tf_from_joint_to_com);

  /**
   * @fn Link getLink(const int id)
   * @brief Get link by id, which is the same as the order of the link added
   *
   * @param id
   * @return Link
   * @sa
   * @detail Detailed description
   */
  const Link &getLink(const int id) const;
  const Link &getBaseLink() const;
  const Link &getEndEffector(const int end_effector_id) const;

  /**
   * @fn Joint getJoint(const int id)
   * @brief Get joint by id, which is the same as the order of the joint added.
   * If the id is -1, return the base_reference_.
   */
  const Joint &getJoint(const int id) const;
  const std::vector<int> &getEndEffectorIdArray() const;
  const std::vector<int> &getLinkIdChainToEndEffector(const int end_effector_id) const;

  /**
   * @fn int getLinkNumber() const
   * @brief Get the number of links in the linkage
   *
   * @return int
   * @sa
   * @detail Detailed description
   */
  int getLinkNumber() const;
  int getJointNumber() const;
  int getEndEffectorNumber() const;
  double getTotalMass() const;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_LINKAGE_HPP_
