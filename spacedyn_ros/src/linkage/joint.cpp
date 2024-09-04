#include "spacedyn_ros/linkage/joint.hpp"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"
#include "spacedyn_ros/linkage/joint_state.hpp"

#include "iostream"

namespace spacedyn_ros {

Joint::Joint(const std::string name, Type type) {
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
  this->tf_to_parent_link_com_ =
      Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
  this->tf_from_parent_link_com_ =
      Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero());
}

void Joint::checkTypeToInput(const Type type) const {
  switch (type) {
  case Type::kRevolute:
    break;
  case Type::kPrismatic:
    break;
  case Type::kUndefined:
    throw std::runtime_error("Error: Joint type is undefined");
  default:
    throw std::invalid_argument("Error: Joint type is incorrect. Please select "
                                "from kFloatingBase. kRevolute, kPrismatic");
  }
}

// Getters
std::string Joint::getName() const { return name_; }
Joint::Type Joint::getType() const { return type_; }
int Joint::getId() const { return id_; }

// Connect joint
void Joint::connect(const int id, const Transform &tf_from_parent_link_com) {
  try { // Check if joint is already connected
    if (id_ != kUndefined) {
      throw std::runtime_error("Error: Joint is already "
                               "connected. Joint id=" +
                               std::to_string(id_));
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
  return;
}

const Transform &Joint::getTransformToParentLinkCom() const { return tf_to_parent_link_com_; }
const Transform &Joint::getTransformFromParentLinkCom() const { return tf_from_parent_link_com_; }

Transform Joint::transformByActuation(const JointState &joint_state) const {
  switch (type_) {
  case Type::kRevolute: {
    Eigen::Matrix3d rot =
        Eigen::AngleAxisd(joint_state.getPosition(), Eigen::Vector3d::UnitZ()).toRotationMatrix();
    return Transform(Frame::kLocal, rot, Eigen::Vector3d::Zero());
  } break;

  case Type::kPrismatic: {
    Eigen::Vector3d trs = Eigen::Vector3d::UnitZ() * joint_state.getPosition();
    return Transform(Frame::kLocal, Eigen::Matrix3d::Identity(), trs);
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
    Eigen::VectorXd twist(6);
    twist.block(0, 0, 3, 1) = Eigen::Vector3d::Zero();
    twist.block(3, 0, 3, 1) = joint_state.getAxisInWorldFrame() * joint_state.getVelocity();
    return Twist(Frame::kWorld, twist);
  } break;

  case Type::kPrismatic: {
    Eigen::VectorXd twist(6);
    twist.block(0, 0, 3, 1) = joint_state.getAxisInWorldFrame() * joint_state.getVelocity();
    twist.block(3, 0, 3, 1) = Eigen::Vector3d::Zero();
    return Twist(Frame::kWorld, twist);
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
    Eigen::VectorXd accel(6);
    accel.block(0, 0, 3, 1) = Eigen::Vector3d::Zero();
    accel.block(3, 0, 3, 1) =
        joint_state.getAxisInWorldFrame() * joint_state.getAcceleration() +
        joint_state.getAxisDerivativeInWorldFrame() * joint_state.getVelocity();
    return Accel(Frame::kWorld, accel);
  } break;

  case Type::kPrismatic: {
    Eigen::VectorXd accel(6);
    accel.block(0, 0, 3, 1) =
        joint_state.getAxisInWorldFrame() * joint_state.getAcceleration() +
        joint_state.getAxisDerivativeInWorldFrame() * joint_state.getVelocity() * 2;
    accel.block(3, 0, 3, 1) = Eigen::Vector3d::Zero();
    return Accel(Frame::kWorld, accel);
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
    return joint_state.getAxisInWorldFrame().dot(
        joint_state.getWrenchToChildInWorldFrame().getOriginTorque());
  } break;

  case Type::kPrismatic: {
    return joint_state.getAxisInWorldFrame().dot(
        joint_state.getWrenchToChildInWorldFrame().getOriginForce());
  } break;

  case Type::kUndefined:
    throw std::invalid_argument("Error: Joint type is undefined");
    break;

  default:
    throw std::logic_error("Unexpected error in Joint::computeEffortToAchieveWrench()");
    break;
  }
}

Eigen::VectorXd Joint::computeVelocityContributionToLinkAccel(const JointState &joint_state,
                                                              const LinkState &link_state) const {
  Eigen::VectorXd jacob_col = Eigen::VectorXd::Zero(6);
  Eigen::Vector3d trans_joint_to_link =
      joint_state.getPoseInWorldFrame().computeTranslationToPoint(link_state.getPoseInWorldFrame());
  Eigen::Vector3d relative_velocity =
      (link_state.getTwistInWorldFrame() - joint_state.getTwistInWorldFrame())
          .getOriginLinierVelocity();

  switch (type_) {
  case Type::kRevolute:
    jacob_col.head(3) = joint_state.getAxisDerivativeInWorldFrame().cross(trans_joint_to_link) +
                        joint_state.getAxisInWorldFrame().cross(relative_velocity);
    jacob_col.tail(3) = joint_state.getAxisDerivativeInWorldFrame();
    return jacob_col;
    break;

  case Joint::Type::kPrismatic:
    jacob_col.head(3) = joint_state.getAxisDerivativeInWorldFrame();
    return jacob_col;
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
