#include "spacedyn_ros/linkage/joint.hpp"
#include "spacedyn_ros/linkage/joint_state.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <iostream>

namespace spacedyn_ros {

Joint::Joint(const std::string name, Type type, const Eigen::Vector3d &axis) {
  try {
    checkTypeToInput(type);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to initialize Joint. Joint type is incorrect.");
  }

  // Initialize
  this->name_ = name;
  this->type_ = type;
  this->id_ = kUndefined;
  this->actuator_id_ = kUndefined;
  this->tf_to_parent_link_com_ =
      Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
  this->tf_from_parent_link_com_ =
      Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
  this->axis_ = axis.normalized();
}

void Joint::checkTypeToInput(const Type type) const {
  switch (type) {
  case Type::kRevolute:
    break;
  case Type::kPrismatic:
    break;
  case Type::kUndefined:
    break;
  case Type::kFixed:
    break;
  default:
    throw std::invalid_argument("Error: Joint type is incorrect.");
  }
}

// Getters
std::string Joint::getName() const { return name_; }
Joint::Type Joint::getType() const { return type_; }
int Joint::getId() const { return id_; }
int Joint::getActuatorId() const { return actuator_id_; }
bool Joint::isActuator() const {
  bool is_actuator = (type_ == Joint::Type::kPrismatic) || (type_ == Joint::Type::kRevolute);
  return is_actuator;
}

// Connect joint
void Joint::connect(const int id, const int actuator_id, const Transform &tf_from_parent_link_com) {
  try { // Check if joint is already connected
    if (id_ != ID::kUndefined) {
      throw std::runtime_error("Error: Joint is already "
                               "connected. Joint id=" +
                               std::to_string(id_));
    }
    if (actuator_id > id) {
      throw std::runtime_error("Error: Actuator id should be lower than joint id.");
    }
    if (type_ == Type::kUndefined) {
      throw std::runtime_error("Error: Joint type is undefined.");
    }
    if (tf_from_parent_link_com.getFrame() != Frame::kLocal) {
      throw std::runtime_error("Error: Transform only in local frame is supported.");
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to connect joint.");
  }
  this->id_ = id;
  this->tf_from_parent_link_com_ = tf_from_parent_link_com;
  this->tf_to_parent_link_com_ = tf_from_parent_link_com.inverse();
  if (isActuator()) {
    this->actuator_id_ = actuator_id;
  }
  return;
}

const Transform &Joint::getTransformToParentLinkCom() const { return tf_to_parent_link_com_; }
const Transform &Joint::getTransformFromParentLinkCom() const { return tf_from_parent_link_com_; }

const Eigen::Vector3d &Joint::getAxisInLocalFrame() const { return axis_; }

Eigen::Vector3d Joint::getAxisInWorldFrame(const JointState &joint_state) const {
  return joint_state.getPoseInWorldFrame().getAttitudeInWorldFrame() * getAxisInLocalFrame();
}

Eigen::Vector3d Joint::getAxisDerivativeInWorldFrame(const JointState &joint_state) const {
  auto angular_velocity = joint_state.getTwistInWorldFrame().getAngularVelocity();
  return angular_velocity.cross(getAxisInWorldFrame(joint_state));
}

Transform Joint::transformByActuation(const JointState &joint_state) const {
  switch (type_) {
  case Type::kRevolute: {
    Eigen::Matrix3d rot =
        Eigen::AngleAxisd(joint_state.getPosition(), getAxisInLocalFrame()).toRotationMatrix();
    return Transform(Frame::kLocal, rot, Eigen::Vector3d::Zero());
  } break;

  case Type::kPrismatic: {
    Eigen::Vector3d trs = getAxisInLocalFrame() * joint_state.getPosition();
    return Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), trs);
  } break;

  case Type::kFixed: {
    return Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
  } break;

  case Type::kUndefined:
    throw std::invalid_argument("Error: Joint type is undefined");
    break;

  default:
    throw std::logic_error("Unexpected error in Joint::computeTransform()");
    break;
  }
}

Twist Joint::twistByActuation(const JointState &joint_state) const {
  switch (type_) {
  case Type::kRevolute: {
    Eigen::Vector6d twist;
    twist.block(0, 0, 3, 1) = Eigen::Vector3d::Zero();
    twist.block(3, 0, 3, 1) = getAxisInWorldFrame(joint_state) * joint_state.getVelocity();
    return Twist(Frame::kWorld, twist);
  } break;

  case Type::kPrismatic: {
    Eigen::Vector6d twist;
    twist.block(0, 0, 3, 1) = getAxisInWorldFrame(joint_state) * joint_state.getVelocity();
    twist.block(3, 0, 3, 1) = Eigen::Vector3d::Zero();
    return Twist(Frame::kWorld, twist);
  } break;

  case Type::kFixed: {
    return Twist(Frame::kWorld, Eigen::Vector6d::Zero());
  } break;

  case Type::kUndefined:
    throw std::invalid_argument("Error: Joint type is undefined");
    break;

  default:
    throw std::logic_error("Unexpected error in Joint::computeTransform()");
    break;
  }
}

Accel Joint::accelByActuation(const JointState &joint_state) const {
  switch (type_) {
  case Type::kRevolute: {
    Eigen::Vector6d accel;
    accel.block(0, 0, 3, 1) = Eigen::Vector3d::Zero();
    accel.block(3, 0, 3, 1) =
        getAxisInWorldFrame(joint_state) * joint_state.getAcceleration() +
        getAxisDerivativeInWorldFrame(joint_state) * joint_state.getVelocity();
    return Accel(Frame::kWorld, accel);
  } break;

  case Type::kPrismatic: {
    Eigen::Vector6d accel;
    accel.block(0, 0, 3, 1) =
        getAxisInWorldFrame(joint_state) * joint_state.getAcceleration() +
        getAxisDerivativeInWorldFrame(joint_state) * joint_state.getVelocity() * 2;
    accel.block(3, 0, 3, 1) = Eigen::Vector3d::Zero();
    return Accel(Frame::kWorld, accel);
  } break;

  case Type::kFixed: {
    return Accel(Frame::kWorld, Eigen::Vector6d::Zero());
  } break;

  case Type::kUndefined:
    throw std::invalid_argument("Error: Joint type is undefined");
    break;

  default:
    throw std::logic_error("Unexpected error in Joint::computeTransform()");
    break;
  }
}

double Joint::computeEffortToAchieveWrench(const JointState &joint_state) const {
  switch (type_) {
  case Type::kRevolute: {
    return getAxisInWorldFrame(joint_state)
        .dot(joint_state.getWrenchToChildInWorldFrame().getTorque());
  } break;

  case Type::kPrismatic: {
    return getAxisInWorldFrame(joint_state)
        .dot(joint_state.getWrenchToChildInWorldFrame().getForce());
  } break;

  case Type::kFixed: {
    return 0.0;
  } break;

  case Type::kUndefined:
    throw std::invalid_argument("Error: Joint type is undefined");
    break;

  default:
    throw std::logic_error("Unexpected error in Joint::computeEffortToAchieveWrench()");
    break;
  }
}

Eigen::Vector6d Joint::computeVelocityContributionToLinkAccel(const JointState &joint_state,
                                                              const LinkState &link_state) const {
  Eigen::Vector6d jacob_col = Eigen::Vector6d::Zero();
  Eigen::Vector3d trans_joint_to_link =
      joint_state.getPoseInWorldFrame().computeTranslationToPoint(link_state.getPoseInWorldFrame());
  Eigen::Vector3d relative_velocity =
      (link_state.getTwistInWorldFrame() - joint_state.getTwistInWorldFrame()).getLinearVelocity();

  switch (type_) {
  case Type::kRevolute:
    jacob_col.head(3) = getAxisDerivativeInWorldFrame(joint_state).cross(trans_joint_to_link) +
                        getAxisInWorldFrame(joint_state).cross(relative_velocity);
    jacob_col.tail(3) = getAxisDerivativeInWorldFrame(joint_state);
    return jacob_col;
    break;

  case Joint::Type::kPrismatic:
    jacob_col.head(3) = getAxisDerivativeInWorldFrame(joint_state);
    return jacob_col;
    break;

  case Joint::Type::kFixed:
    return Eigen::Vector6d::Zero(6);
    break;

  case Type::kUndefined:
    throw std::invalid_argument("Error: Joint type is undefined");
    break;

  default:
    throw std::logic_error("Unexpected error in Joint::computeRelativeTwistEffect()");
    break;
  }
}
} // namespace spacedyn_ros
