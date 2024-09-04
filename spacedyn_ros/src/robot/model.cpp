#include "spacedyn_ros/robot/model.hpp"
#include "eigen3/Eigen/Core"
#include "urdf/model.h"

#include "iostream"

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

  std::map<std::string, int> link_name_to_index;
  int link_index = 0;

  std::vector<Link> spd_links;

  // Load links with mass and inertia and map link names to indices
  for (const auto &link_pair : urdf_model.links_) {
    const auto &link = link_pair.second;

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

    spd_links.push_back(Link(link->name, Inertia(Frame::kLocal, mass, inertia_matrix)));
    link_name_to_index[link->name] = link_index++;
  }

  // Load the base link
  linkage_.addBase(spd_links.at(link_name_to_index.at(urdf_model.getRoot()->name)));

  // Load joints and link them with their parent and child links
  for (const auto &joint_pair : urdf_model.joints_) {
    const auto &joint = joint_pair.second;

    Joint::Type joint_type;
    if (joint->type == urdf::Joint::REVOLUTE) {
      joint_type = Joint::Type::kRevolute;
    } else if (joint->type == urdf::Joint::PRISMATIC) {
      joint_type = Joint::Type::kPrismatic;
    } else {
      throw std::runtime_error("Unsupported joint type.");
    }

    Joint spd_joint(joint->name, joint_type);

    // Find the parent and child link indices
    int parent_link_index = link_name_to_index[joint->parent_link_name];
    int child_link_index = link_name_to_index[joint->child_link_name];

    // Set the joint transforms
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

    Transform tf_parent_link_to_joint(Frame::kLocal, parent_joint_to_parent_link.inverse() *
                                                         parent_joint_to_joint);
    Transform tf_joint_to_child_link(Frame::kLocal, joint_to_child_link);

    // Add the joint and link together
    linkage_.addJointWithLink(parent_link_index, spd_joint, tf_parent_link_to_joint,
                              spd_links.at(child_link_index), // Assume child link already created
                              tf_joint_to_child_link);
  }
}

const Linkage &Model::getLinkage() const { return this->linkage_; }
double Model::getTotalMass() const { return getLinkage().getTotalMass(); }
int Model::getLinkNumber() const { return getLinkage().getLinkNumber(); }
int Model::getJointNumber() const { return getLinkage().getJointNumber(); }
double Model::getDeltaTimeSec() const { return this->dt_microsec_ / 1e6; }
double Model::getDeltaTimeMilliSec() const { return this->dt_microsec_ / 1e3; }
size_t Model::getDeltaTimeMicroSec() const { return this->dt_microsec_; }
const Eigen::Vector3d &Model::getGravity() const { return this->gravity_; }

void Model::setGravity(const Eigen::Vector3d &gravity) { this->gravity_ = gravity; }
void Model::setDeltaTimeMicroSec(const size_t dt_microsec) { this->dt_microsec_ = dt_microsec; }
} // namespace spacedyn_ros
