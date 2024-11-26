#include "spacedyn_ros/robot/model.hpp"
#include "urdf/model.h"
#include <eigen3/Eigen/Core>
#include <iostream>

namespace spacedyn_ros {
Model::Model(const Linkage &linkage) : linkage_(linkage) {
  // Default values
  this->dt_microsec_ = 1000;
  this->gravity_ = Eigen::Vector3d(0, 0, -9.8);
}

Model::Model(const std::string &path_to_urdf) {
  // Default values
  // TODO: Load URDF file and create Linkage object
  loadURDF(path_to_urdf);

  // Set other properties
  setGravity(Eigen::Vector3d(0, 0, -9.8)); // default gravity vector
  setDeltaTimeMicroSec(1000);              // default time step
}

// FIXME: Implement this function. Joint axis is not considered yet.
// TODO: Update joint class to include joint axis
void Model::loadURDF(const std::string &path_to_urdf) {
  urdf::Model urdf_model;
  if (!urdf_model.initFile(path_to_urdf)) {
    throw std::runtime_error("Failed to parse URDF file.");
  }

  // ------------------------ Arrange links in the order of the tree ------------------------ //
  // Create a map to store the link index and name
  // Stored from end link to base link
  std::map<int, std::string> link_index_to_name;
  std::map<std::string, int> link_name_to_index;
  std::map<std::string, int> joint_name_to_index;

  auto link = urdf_model.getRoot();
  int link_index = urdf_model.links_.size() - 1; // link index equals to the link_id in spacedyn

  // Variables to traverse the tree
  std::vector<std::string> link_to_be_checked; // name of the link to be checked
  std::vector<int> children_checked_counts;    // children_checked_counts.at(i) is the num of
                                               // children checked for the link_to_be_checked.at(i)

  // Start from the root link
  children_checked_counts.push_back(0);
  link_to_be_checked.push_back(link->name);
  while (true) {
    if (link_to_be_checked.size() == 0) {
      // End of link list. All links have been checked
      break;
    }

    // link to be checked its child links
    link = urdf_model.links_.at(link_to_be_checked.back());

    int uncheck_child_count = link->child_links.size() - children_checked_counts.back();
    if (uncheck_child_count <= 0) {
      // No more child links to check
      // Remove the link from the list and go back to the parent link
      link_index_to_name[link_index] = link->name;
      link_name_to_index[link->name] = link_index;
      if (link_index > 0) {
        // If the link is not the base link, store the parent joint index
        joint_name_to_index[link->parent_joint->name] = link_index - 1;
      }
      link_index--;
      children_checked_counts.pop_back();
      link_to_be_checked.pop_back();

      // Increment link_count means that the link has been checked
      if (!children_checked_counts.empty()) {
        children_checked_counts.back() = children_checked_counts.back() + 1;
      }
      continue;
    }

    // If there are child links to check, add the child link to the list
    link_to_be_checked.push_back(link->child_links.at(children_checked_counts.back())->name);
    children_checked_counts.push_back(0);
  }

  // ------------------------ Create SpaceDyn linkage ------ ------------------------ //
  for (size_t li = 0; li < urdf_model.links_.size(); li++) {
    // ----- Create Link ----- //
    const auto &link = urdf_model.links_.at(link_index_to_name[li]);

    double mass = 0.0;
    Eigen::Matrix3d inertia_matrix = Eigen::Matrix3d::Zero();

    if (link->inertial) {
      mass = link->inertial->mass;

      inertia_matrix(0, 0) = link->inertial->ixx;
      inertia_matrix(1, 1) = link->inertial->iyy;
      inertia_matrix(2, 2) = link->inertial->izz;

      inertia_matrix(0, 1) = inertia_matrix(1, 0) = link->inertial->ixy;
      inertia_matrix(0, 2) = inertia_matrix(2, 0) = link->inertial->ixz;
      inertia_matrix(1, 2) = inertia_matrix(2, 1) = link->inertial->iyz;
    }

    auto spd_link = Link(link->name, Inertia(Frame::kLocal, mass, inertia_matrix));
    if (li == 0) {
      // Base does not have a parent joint
      linkage_.addBase(spd_link);
      continue;
    }

    // ----- Create Joint ----- //
    // TODO: Remove this part when joint axis is implemented
    const auto &joint = urdf_model.joints_.at(link->parent_joint->name);
    bool axis_is_valid = (joint->axis.x == 0.0 && joint->axis.y == 0.0 && joint->axis.z == 1.0);
    if (!axis_is_valid) {
      throw std::runtime_error(
          "Unsupported joint axis. Currently only z-axis is supported for joint axis: " +
          joint->name);
    }
    Joint::Type joint_type;
    if (joint->type == urdf::Joint::REVOLUTE) {
      joint_type = Joint::Type::kRevolute;
    } else if (joint->type == urdf::Joint::PRISMATIC) {
      joint_type = Joint::Type::kPrismatic;
    } else {
      throw std::runtime_error("Unsupported joint type: " + joint->name);
    }

    auto spd_joint = Joint(joint->name, joint_type);

    // ----- Set Joint and Link ----- //
    int parent_link_index = link_name_to_index[joint->parent_link_name];

    // Convert URDF transforms to Eigen::Isometry3d for SpaceDyn
    //
    // URDF:parent_joint -> parent_link
    //      parent_joint -> joint -> child_link
    // => SpaceDyn: parent_link -> joint -> child_link
    Eigen::Vector3d parent_joint_to_joint_trs(joint->parent_to_joint_origin_transform.position.x,
                                              joint->parent_to_joint_origin_transform.position.y,
                                              joint->parent_to_joint_origin_transform.position.z);
    Eigen::Quaterniond parent_joint_to_joint_rot(
        joint->parent_to_joint_origin_transform.rotation.w,
        joint->parent_to_joint_origin_transform.rotation.x,
        joint->parent_to_joint_origin_transform.rotation.y,
        joint->parent_to_joint_origin_transform.rotation.z);
    Eigen::Isometry3d parent_joint_to_joint =
        Eigen::Translation3d(parent_joint_to_joint_trs) * parent_joint_to_joint_rot;

    Eigen::Vector3d parent_joint_to_parent_link_trs(
        urdf_model.links_.at(joint->parent_link_name)->inertial->origin.position.x,
        urdf_model.links_.at(joint->parent_link_name)->inertial->origin.position.y,
        urdf_model.links_.at(joint->parent_link_name)->inertial->origin.position.z);
    Eigen::Quaterniond parent_joint_to_parent_link_rot(
        urdf_model.links_.at(joint->parent_link_name)->inertial->origin.rotation.w,
        urdf_model.links_.at(joint->parent_link_name)->inertial->origin.rotation.x,
        urdf_model.links_.at(joint->parent_link_name)->inertial->origin.rotation.y,
        urdf_model.links_.at(joint->parent_link_name)->inertial->origin.rotation.z);
    Eigen::Isometry3d parent_joint_to_parent_link =
        Eigen::Translation3d(parent_joint_to_parent_link_trs) * parent_joint_to_parent_link_rot;

    Eigen::Vector3d joint_to_child_link_trs(
        urdf_model.links_.at(joint->child_link_name)->inertial->origin.position.x,
        urdf_model.links_.at(joint->child_link_name)->inertial->origin.position.y,
        urdf_model.links_.at(joint->child_link_name)->inertial->origin.position.z);
    Eigen::Quaterniond joint_to_child_link_rot(
        urdf_model.links_.at(joint->child_link_name)->inertial->origin.rotation.w,
        urdf_model.links_.at(joint->child_link_name)->inertial->origin.rotation.x,
        urdf_model.links_.at(joint->child_link_name)->inertial->origin.rotation.y,
        urdf_model.links_.at(joint->child_link_name)->inertial->origin.rotation.z);
    Eigen::Isometry3d joint_to_child_link =
        Eigen::Translation3d(joint_to_child_link_trs) * joint_to_child_link_rot;

    // Set the joint and link transforms
    Transform tf_parent_link_to_joint(Frame::kLocal, parent_joint_to_parent_link.inverse() *
                                                         parent_joint_to_joint);
    Transform tf_joint_to_child_link(Frame::kLocal, joint_to_child_link);

    // Add the joint and link together
    linkage_.addJointWithLink(parent_link_index, spd_joint, tf_parent_link_to_joint, spd_link,
                              tf_joint_to_child_link);
  }
}

const Linkage &Model::getLinkage() const { return this->linkage_; }
double Model::getTotalMass() const { return getLinkage().getTotalMass(); }

int Model::getDof() const { return getLinkage().getDof(); }
int Model::getLinkNumber() const { return getLinkage().getLinkNumber(); }
int Model::getJointNumber() const { return getLinkage().getJointNumber(); }
int Model::getActuatorNumber() const { return getLinkage().getActuatorNumber(); }
int Model::getEndEffectorNumber() const { return getLinkage().getEndEffectorNumber(); }
double Model::getDeltaTimeSec() const { return this->dt_microsec_ / 1e6; }
double Model::getDeltaTimeMilliSec() const { return this->dt_microsec_ / 1e3; }
size_t Model::getDeltaTimeMicroSec() const { return this->dt_microsec_; }
const Eigen::Vector3d &Model::getGravity() const { return this->gravity_; }

void Model::setGravity(const Eigen::Vector3d &gravity) { this->gravity_ = gravity; }
void Model::setDeltaTimeMicroSec(const size_t dt_microsec) { this->dt_microsec_ = dt_microsec; }
} // namespace spacedyn_ros
