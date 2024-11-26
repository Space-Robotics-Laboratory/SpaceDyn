#include "spacedyn_ros/robot/state_variable.hpp"
#include "spacedyn_ros/motion/kinematics.hpp"
#include "spacedyn_ros/robot/model.hpp"
#include <eigen3/Eigen/Core>
#include <iostream>

namespace spacedyn_ros {

StateVariable::StateVariable(Model const &model) {
  try {
    checkModelToSet(model);
    this->link_number_ = model.getLinkage().getLinkNumber();
    joint_number_ = model.getLinkage().getJointNumber();
    actuator_number_ = model.getActuatorNumber();

    this->links_state_.resize(this->link_number_);
    this->joints_state_.resize(joint_number_);
    this->joint_is_actuator_.resize(joint_number_);

    // Init joint and link state
    if (this->link_number_ == 0 && this->joint_number_ == 0) {
      return;
    }
    this->links_state_.at(0) = LinkState();
    for (int link_id = 1; link_id < this->link_number_; link_id++) {
      int joint_id = link_id - 1;
      this->links_state_.at(link_id) = LinkState();
      this->joints_state_.at(joint_id) = JointState();
      this->joint_is_actuator_.at(joint_id) = model.getLinkage().getJoint(joint_id).isActuator();
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to initialize state variable.");
  }
}

StateVariable::StateVariable(const int link_number, const std::vector<bool> &joint_is_actuator) {
  try {
    if (link_number < 1) {
      throw std::invalid_argument("Error: Link number is invalid");
    }
    this->link_number_ = link_number;
    joint_number_ = link_number - 1;
    this->joint_is_actuator_ = joint_is_actuator;

    this->links_state_.resize(this->link_number_);
    this->joints_state_.resize(joint_number_);

    // Init joint and link state
    this->actuator_number_ = 0;
    this->links_state_.at(0) = LinkState();
    for (int link_id = 1; link_id < this->link_number_; link_id++) {
      int joint_id = link_id - 1;
      this->links_state_.at(link_id) = LinkState();
      this->joints_state_.at(joint_id) = JointState();
      if (joint_is_actuator.at(joint_id)) {
        this->actuator_number_++;
      }
    }
    if (actuator_number_ > joint_number_) {
      throw std::invalid_argument("Error: Actuator number is invalid");
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to initialize state variable.");
  }
}

void StateVariable::checkModelToSet(const Model &model) const {
  const int link_number = model.getLinkage().getLinkNumber();
  const int joint_number = model.getLinkage().getJointNumber();
  if (link_number == 0 && joint_number == 0) {
    // HACK: Allow empty model to create robot before setting model
    return;
  } else if (link_number != joint_number + 1) {
    throw std::logic_error("Error:Model is invalid. Joint number=" + std::to_string(joint_number) +
                           "is not equal to link number - 1=" + std::to_string(link_number - 1));
  }
}

void StateVariable::checkLinkIdToCall(const int link_id) const {
  if (link_id < 0 || link_id >= getLinkNumber()) {
    throw std::out_of_range("Link ID out of range");
  }
}

void StateVariable::checkJointIdToCall(const int joint_id) const {
  if (joint_id < 0 || joint_id >= getJointNumber()) {
    throw std::out_of_range("Joint ID out of range");
  }
}

void StateVariable::checkJointVectorToSet(const Eigen::VectorXd &joint_vector) const {
  if (joint_vector.size() != actuator_number_) {
    throw std::invalid_argument("Error: Joint vector size is invalid");
  }
}

void StateVariable::checkGeneralizedVectorToSet(const Eigen::VectorXd &generalized_vector) const {
  if (generalized_vector.size() != getDof()) {
    throw std::invalid_argument("Error: Generalized vector size is invalid");
  }
}

StateVariable StateVariable::copyVacant() const {
  StateVariable state_variable(this->link_number_, this->joint_is_actuator_);
  return state_variable;
}

int StateVariable::getLinkNumber() const { return this->link_number_; }
int StateVariable::getJointNumber() const { return joint_number_; }
int StateVariable::getActuatorNumber() const { return actuator_number_; }

int StateVariable::getDof() const { return actuator_number_ + 6; }

const LinkState &StateVariable::getLinkState(const int link_id) const {
  try {
    checkLinkIdToCall(link_id);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to get link state.");
  }
  return this->links_state_.at(link_id);
}

const LinkState &StateVariable::getLinkState(const Link &link) const {
  return this->links_state_.at(link.getId());
}

const JointState &StateVariable::getJointState(const int joint_id) const {
  try {
    checkJointIdToCall(joint_id);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to get joint state.");
  }

  return this->joints_state_.at(joint_id);
}

const JointState &StateVariable::getJointState(const Joint &joint) const {
  return this->joints_state_.at(joint.getId());
}

Eigen::VectorXd StateVariable::getJointPosition() const {
  Eigen::VectorXd joint_position(actuator_number_);
  int actuator_id = 0;
  for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
    if (this->joint_is_actuator_.at(joint_id)) {
      joint_position(actuator_id) = this->joints_state_.at(joint_id).getPosition();
      actuator_id++;
    }
  }
  return joint_position;
}

Eigen::VectorXd StateVariable::getJointVelocity() const {
  Eigen::VectorXd joint_velocity(actuator_number_);
  int actuator_id = 0;
  for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
    if (this->joint_is_actuator_.at(joint_id)) {
      joint_velocity(actuator_id) = this->joints_state_.at(joint_id).getVelocity();
      actuator_id++;
    }
  }
  return joint_velocity;
}

Eigen::VectorXd StateVariable::getJointAcceleration() const {
  Eigen::VectorXd joint_acceleration(actuator_number_);
  int actuator_id = 0;
  for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
    if (this->joint_is_actuator_.at(joint_id)) {
      joint_acceleration(actuator_id) = this->joints_state_.at(joint_id).getAcceleration();
      actuator_id++;
    }
  }
  return joint_acceleration;
}

Eigen::VectorXd StateVariable::getJointEffort() const {
  Eigen::VectorXd joint_effort(actuator_number_);
  int actuator_id = 0;
  for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
    if (this->joint_is_actuator_.at(joint_id)) {
      joint_effort(actuator_id) = this->joints_state_.at(joint_id).getEffort();
      actuator_id++;
    }
  }
  return joint_effort;
}

Eigen::VectorXd StateVariable::getGeneralizedVelocity() const {
  Eigen::VectorXd generalized_velocity(getDof());
  generalized_velocity.head(6) =
      this->links_state_.at(Link::ID::kBase).getTwistInWorldFrame().getTwist();
  generalized_velocity.tail(actuator_number_) = getJointVelocity();
  return generalized_velocity;
}

Eigen::VectorXd StateVariable::getGeneralizedAcceleration() const {
  Eigen::VectorXd generalized_acceleration(getDof());
  generalized_acceleration.head(6) =
      this->links_state_.at(Link::ID::kBase).getAccelInWorldFrame().getAccel();
  generalized_acceleration.tail(actuator_number_) = getJointAcceleration();
  return generalized_acceleration;
}

void StateVariable::setLinkPoseInWorldFrame(const int link_id, const Eigen::Isometry3d &pose) {
  try {
    checkLinkIdToCall(link_id);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to set link pose.");
  }
  this->links_state_.at(link_id).setPoseInWorldFrame(Pose(pose));
}

void StateVariable::setLinkTwistInWorldFrame(const int link_id, const Eigen::Vector6d &twist) {
  try {
    checkLinkIdToCall(link_id);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to set link twist.");
  }
  this->links_state_.at(link_id).setTwistInWorldFrame(Twist(Frame::kWorld, twist));
}

void StateVariable::setLinkAccelInWorldFrame(const int link_id, const Eigen::Vector6d &accel) {
  try {
    checkLinkIdToCall(link_id);
    this->links_state_.at(link_id).setAccelInWorldFrame(Accel(Frame::kWorld, accel));
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to set link accel.");
  }
}

void StateVariable::setLinkExternallyAppliedWrenchInWorldFrame(const int link_id,
                                                               const Eigen::VectorXd &wrench) {
  try {
    checkLinkIdToCall(link_id);
    this->links_state_.at(link_id).setExternallyAppliedWrenchInWorldFrame(
        Wrench(Frame::kWorld, wrench));
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to set link external wrench.");
  }
}

void StateVariable::setJointPosition(const Eigen::VectorXd &joint_position) {
  try {
    checkJointVectorToSet(joint_position);
    int actuator_id = 0;
    for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
      if (this->joint_is_actuator_.at(joint_id)) {
        this->joints_state_.at(joint_id).setPosition(joint_position(actuator_id));
        actuator_id++;
      }
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set joint position.");
  }
}

void StateVariable::setJointVelocity(const Eigen::VectorXd &joint_velocity) {
  try {
    checkJointVectorToSet(joint_velocity);
    int actuator_id = 0;
    for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
      if (this->joint_is_actuator_.at(joint_id)) {
        this->joints_state_.at(joint_id).setVelocity(joint_velocity(actuator_id));
        actuator_id++;
      }
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set joint velocity.");
  }
}

void StateVariable::setJointAcceleration(const Eigen::VectorXd &joint_acceleration) {
  try {
    checkJointVectorToSet(joint_acceleration);
    int actuator_id = 0;
    for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
      if (this->joint_is_actuator_.at(joint_id)) {
        this->joints_state_.at(joint_id).setAcceleration(joint_acceleration(actuator_id));
        actuator_id++;
      }
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set joint acceleration.");
  }
}

void StateVariable::setJointEffort(const Eigen::VectorXd &joint_effort) {
  try {
    checkJointVectorToSet(joint_effort);
    int actuator_id = 0;
    for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
      if (this->joint_is_actuator_.at(joint_id)) {
        this->joints_state_.at(joint_id).setEffort(joint_effort(actuator_id));
        actuator_id++;
      }
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set joint effort.");
  }
}

void StateVariable::setLinkState(const int link_id, const LinkState &link_state) {
  try {
    checkLinkIdToCall(link_id);
    this->links_state_.at(link_id) = link_state;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to set link state.");
  }
}

void StateVariable::setJointState(const int joint_id, const JointState &joint_state) {
  try {
    checkJointIdToCall(joint_id);
    if (!this->joint_is_actuator_.at(joint_id) && joint_state.hasNonZeroState()) {
      throw std::invalid_argument("Error: Joint is not actuator but has non-zero state.");
    }
    this->joints_state_.at(joint_id) = joint_state;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to set joint state.");
  }
}

void StateVariable::setGeneralizedCoordinates(const Eigen::Isometry3d &base_pose,
                                              const Eigen::VectorXd &joint_position) {
  try {
    setLinkPoseInWorldFrame(Link::ID::kBase, base_pose);
    setJointPosition(joint_position);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set generalized coordinates.");
  }
}

void StateVariable::setGeneralizedVelocity(const Eigen::VectorXd &generalized_velocity) {
  const int DOF = 6;
  try {
    checkGeneralizedVectorToSet(generalized_velocity);
    setLinkTwistInWorldFrame(Link::ID::kBase, generalized_velocity.head(DOF));
    setJointVelocity(generalized_velocity.tail(actuator_number_));
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set generalized velocity.");
  }
}

void StateVariable::setGeneralizedAcceleration(const Eigen::VectorXd &generalized_acceleration) {
  const int DOF = 6;
  try {
    checkGeneralizedVectorToSet(generalized_acceleration);
    setLinkAccelInWorldFrame(Link::ID::kBase, generalized_acceleration.head(DOF));
    setJointAcceleration(generalized_acceleration.tail(actuator_number_));
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set generalized acceleration.");
  }
}
} // namespace spacedyn_ros
