#include "spacedyn_ros/linkage/link.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <iostream>

namespace spacedyn_ros {
Link::Link(const std::string name, const Inertia &inertia_in_local_frame,
           const Transform &tf_com_to_end_tip) {
  try {
    if (inertia_in_local_frame.getFrame() != Frame::kLocal) {
      throw std::invalid_argument("Error: Inertia must be expressed in local frame.");
    }
    if (tf_com_to_end_tip.getFrame() != Frame::kLocal) {
      throw std::invalid_argument("Error: Transform must be expressed in local frame.");
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to create link.");
  }

  // Assignment
  this->name_ = name;
  this->inertia_in_local_frame_ = inertia_in_local_frame;
  this->tf_from_com_to_parent_joint_ =
      Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
  this->tf_from_parent_joint_to_com_ =
      Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());

  // Initialization
  this->id_ = kUndefined;
  this->parent_id_ = kUndefined;
  this->children_number_ = 0;
  this->is_end_effector_ = false;
  this->tf_com_to_end_tip_ = tf_com_to_end_tip;
}

void Link::connect(const int parent_id, const int assigned_id,
                   const Transform &tf_from_parent_joint_to_com) {
  try {
    if (id_ != kUndefined) {
      throw std::runtime_error("Error: Link is already "
                               "connected. Link id=" +
                               std::to_string(id_));
    }
    if (tf_from_parent_joint_to_com.getFrame() != Frame::kLocal) {
      throw std::runtime_error("Error: Transform only in local frame is supported.");
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to connect link.");
  }

  this->id_ = assigned_id;
  this->parent_id_ = parent_id;
  this->tf_from_parent_joint_to_com_ = tf_from_parent_joint_to_com;
  this->tf_from_com_to_parent_joint_ = tf_from_parent_joint_to_com.inverse();
  this->is_end_effector_ = true;
  return;
}

void Link::acceptChild() {
  this->children_number_++;
  this->is_end_effector_ = false;
  // Transform to end tip has no meaning
  this->tf_com_to_end_tip_ =
      Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
  return;
}

std::string Link::getName() const { return name_; }

int Link::getId() const {
  if (id_ == kUndefined) {
    throw std::runtime_error("Error: Failed to get link id. Link is not connected.");
  }

  return id_;
}

int Link::getParentId() const { return parent_id_; }

int Link::getParentJointId() const {
  if (id_ == kBase) {
    throw std::runtime_error("Error: Failed to get parent joint id. Link is a base link.");
  }

  return id_ - 1;
}

int Link::getChildrenNumber() const { return children_number_; }

double Link::getMass() const { return inertia_in_local_frame_.getMass(); }

const Inertia &Link::getInertiaInLocalFrame() const { return inertia_in_local_frame_; }

Inertia Link::computeInertiaInWorldFrame(const Pose &pose) const {
  return inertia_in_local_frame_.getInertiaInFrame(Frame::kWorld, pose);
}

const Transform &Link::getTransformToParentJoint() const { return tf_from_com_to_parent_joint_; }
const Transform &Link::getTransformFromParentJoint() const { return tf_from_parent_joint_to_com_; }
const Transform &Link::getTransformToEndTip() const {
  if (!is_end_effector_) {
    throw std::runtime_error("Error: Link is not an end effector.");
  }
  return tf_com_to_end_tip_;
}

bool Link::isBase() const { return id_ == kBase; }
bool Link::isEndEffector() const { return is_end_effector_; }
} // namespace spacedyn_ros
