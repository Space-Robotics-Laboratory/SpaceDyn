#include "spacedyn_ros/linkage/linkage.hpp"
#include "spacedyn_ros/linkage/joint.hpp"
#include "spacedyn_ros/linkage/link.hpp"
#include <iostream>

namespace spacedyn_ros {

// TODO: Make Connection class to store parent and child id
Linkage::Linkage() {
  // Initialize
  this->total_mass_ = 0;
  this->joints_.resize(0);
  this->links_.resize(0);
  this->actuators_.resize(0);
  this->end_effectors_.resize(0);
}

void Linkage::checkLinkIdToCall(const int link_id) const {
  if (link_id < 0 || link_id >= getLinkNumber()) {
    throw std::out_of_range("Error: Id not found. Id=" + std::to_string(link_id) +
                            " should be lower than link_number=" + std::to_string(getLinkNumber()) +
                            " and non-negative");
  }
}

void Linkage::checkJointIdToCall(const int joint_id) const {
  if (joint_id < 0 || joint_id >= getJointNumber()) {
    throw std::out_of_range("Error: Id not found. Id=" + std::to_string(joint_id) +
                            " should be lower than joint_number=" +
                            std::to_string(getJointNumber()) + " and non-negative");
  }
}

void Linkage::checkEndEffectorIdToCall(const int end_effector_id) const {
  if (end_effector_id < 0 || end_effector_id >= getEndEffectorNumber()) {
    throw std::out_of_range("Error: Id not found. Id=" + std::to_string(end_effector_id) +
                            " should be lower than end_effector_number=" +
                            std::to_string(getEndEffectorNumber()) + " and positive");
  }
}

void Linkage::addBase(const Link &link) {
  // Check if base is already added
  if (getLinkNumber() != 0) {
    throw std::runtime_error("Error: Failed to connect base. There is already "
                             "a Base link.");
  }

  // Copy link to be added
  auto base_link = link;

  // Set base link
  base_link.connect(Link::kUndefined, Link::kBase,
                    Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero()));

  // Add end effector
  this->end_effectors_.push_back(base_link.getId());

  // Add link to linkage
  links_.push_back(base_link);

  // Update total mass
  this->total_mass_ += base_link.getInertiaInLocalFrame().getMass();

  // Update connectivity to end effector

  // Set link connection graph of base
  connection_.addBase();
}

void Linkage::addLink(const int parent_link_id, const Link &child_input,
                      const Transform &tf_from_parent_joint_to_com) {
  // Copy link to be added
  auto child_link = child_input;

  // Set id
  const int child_link_id_to_set = getLinkNumber();
  child_link.connect(parent_link_id, child_link_id_to_set, tf_from_parent_joint_to_com);
  // Set parent as non-end effector
  auto parent_link = getLink(parent_link_id);

  // Reconstruction of end effector
  this->end_effectors_ = replaceEndEffector(parent_link, child_link);
  parent_link.acceptChild();
  replaceLink(parent_link_id, parent_link);

  // Add link to linkage
  links_.push_back(child_link);

  // Update total mass
  this->total_mass_ += child_link.getInertiaInLocalFrame().getMass();

  // Update connectivity to end effector

  // Set link connection graph of link
  connection_.addLink(parent_link_id, child_link_id_to_set);
}

void Linkage::replaceLink(const int id, const Link &link) {
  // update total mass
  this->total_mass_ -= this->links_.at(id).getInertiaInLocalFrame().getMass();
  this->total_mass_ += link.getInertiaInLocalFrame().getMass();

  // Replace link
  this->links_.at(id) = link;
}

void Linkage::addJoint(const Joint &joint_input,
                       const Transform &tf_from_parent_link_com_to_joint) {
  // Copy joint to be added
  auto child_joint = joint_input;

  // Set id
  child_joint.connect(getJointNumber(), getActuatorNumber(), tf_from_parent_link_com_to_joint);

  // Add link to linkage
  joints_.push_back(child_joint);

  // If it's actuator, add to actuator list
  if (child_joint.isActuator()) {
    actuators_.push_back(child_joint.getId());
  }
}

void Linkage::addJointWithLink(const int parent_link_id, const Joint &child_joint,
                               const Transform &tf_from_parent_link_com, const Link &child_link,
                               const Transform &tf_from_joint_to_com) {
  try {
    checkLinkIdToCall(parent_link_id);
    if (getLinkNumber() == 0) {
      throw std::invalid_argument("Error: Base does not have a joint. Please add "
                                  "Base using addBase() function before add Joint.");
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to connect joint with link.");
  }

  // Add joint
  addJoint(child_joint, tf_from_parent_link_com);

  // Add link
  addLink(parent_link_id, child_link, tf_from_joint_to_com);
}

const Link &Linkage::getLink(const int id) const {
  try {
    checkLinkIdToCall(id);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to get link.");
  }

  return links_.at(id);
}

const Link &Linkage::getBaseLink() const {
  try {
    checkLinkIdToCall(Link::kBase);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to get base link.");
  }

  return links_.at(Link::kBase);
}

const Link &Linkage::getEndEffector(const int end_effector_id) const {
  try {
    checkEndEffectorIdToCall(end_effector_id);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to get end effector.");
  }

  int link_id = end_effectors_.at(end_effector_id);

  return getLink(link_id);
}

const Joint &Linkage::getJoint(const int id) const {
  try {
    checkJointIdToCall(id);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to get joint.");
  }

  return joints_.at(id);
}

const std::vector<int> &Linkage::getEndEffectorIdArray() const { return end_effectors_; }

std::vector<int> Linkage::getLinkIdChain(const int start_link_id, const int end_link_id) const {
  return connection_.getLinkIdChain(start_link_id, end_link_id);
}

std::vector<int> Linkage::replaceEndEffector(const Link parent, const Link child) const {
  std::vector<int> end_effectors;
  end_effectors = this->end_effectors_;

  // If the parent is end effector, replace its id in end_effectors vector
  // with child's id
  if (parent.isEndEffector()) {
    for (int i = 0; i < static_cast<int>(end_effectors_.size()); i++) {
      if (end_effectors.at(i) == parent.getId()) {
        end_effectors.at(i) = child.getId();
      }
    }
  } else {
    end_effectors.push_back(child.getId());
  }

  return end_effectors;
}

int Linkage::getDof() const { return getActuatorNumber() + 6; }
int Linkage::getLinkNumber() const { return links_.size(); }
int Linkage::getJointNumber() const { return joints_.size(); }
int Linkage::getActuatorNumber() const { return actuators_.size(); }
int Linkage::getEndEffectorNumber() const { return end_effectors_.size(); }
double Linkage::getTotalMass() const { return this->total_mass_; }
} // namespace spacedyn_ros
