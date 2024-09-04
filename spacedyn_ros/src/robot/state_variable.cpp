#include "spacedyn_ros/robot/state_variable.hpp"
#include "eigen3/Eigen/Core"
#include "iostream"
#include "spacedyn_ros/motion/kinematics.hpp"
#include "spacedyn_ros/robot/model.hpp"

namespace spacedyn_ros {

StateVariable::StateVariable(Model const &model) {
  try {
    checkModelToSet(model);
    this->link_number_ = model.getLinkage().getLinkNumber();
    joint_number_ = model.getLinkage().getJointNumber();

    this->links_state_.resize(this->link_number_);
    this->joints_state_.resize(joint_number_);

    // Init joint and link state
    if (this->link_number_ == 0 && this->joint_number_ == 0) {
      return;
    }
    this->links_state_.at(0) = LinkState();
    for (int link_id = 1; link_id < this->link_number_; link_id++) {
      int joint_id = link_id - 1;
      this->links_state_.at(link_id) = LinkState();
      this->joints_state_.at(joint_id) = JointState();
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to initialize state variable.");
  }
}

StateVariable::StateVariable(const int link_number) {
  try {
    if (link_number < 1) {
      throw std::invalid_argument("Error: Link number is invalid");
    }
    this->link_number_ = link_number;
    joint_number_ = link_number - 1;

    this->links_state_.resize(this->link_number_);
    this->joints_state_.resize(joint_number_);

    // Init joint and link state
    this->links_state_.at(0) = LinkState();
    for (int link_id = 1; link_id < this->link_number_; link_id++) {
      int joint_id = link_id - 1;
      this->links_state_.at(link_id) = LinkState();
      this->joints_state_.at(joint_id) = JointState();
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
    // Vacant model
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
  if (joint_vector.size() != joint_number_) {
    throw std::invalid_argument("Error: Joint vector size is invalid");
  }
}

void StateVariable::checkGeneralizedVectorToSet(const Eigen::VectorXd &generalized_vector) const {
  const int DOF = 6;
  if (generalized_vector.size() != joint_number_ + DOF) {
    throw std::invalid_argument("Error: Generalized vector size is invalid");
  }
}

StateVariable StateVariable::copyVacant() const {
  StateVariable state_variable(this->link_number_);
  return state_variable;
}

int StateVariable::getLinkNumber() const { return this->link_number_; }
int StateVariable::getJointNumber() const { return joint_number_; }

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
  Eigen::VectorXd joint_position(joint_number_);
  for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
    joint_position(joint_id) = this->joints_state_.at(joint_id).getPosition();
  }
  return joint_position;
}

Eigen::VectorXd StateVariable::getJointVelocity() const {
  Eigen::VectorXd joint_velocity(joint_number_);
  for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
    joint_velocity(joint_id) = this->joints_state_.at(joint_id).getVelocity();
  }
  return joint_velocity;
}

Eigen::VectorXd StateVariable::getJointAcceleration() const {
  Eigen::VectorXd joint_acceleration(joint_number_);
  for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
    joint_acceleration(joint_id) = this->joints_state_.at(joint_id).getAcceleration();
  }
  return joint_acceleration;
}

Eigen::VectorXd StateVariable::getJointEffort() const {
  Eigen::VectorXd joint_effort(joint_number_);
  for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
    joint_effort(joint_id) = this->joints_state_.at(joint_id).getEffort();
  }
  return joint_effort;
}

Eigen::VectorXd StateVariable::getGeneralizedVelocity() const {
  const int DOF = 6;
  Eigen::VectorXd generalized_velocity(joint_number_ + DOF);
  generalized_velocity.head(DOF) =
      this->links_state_.at(Link::ID::kBase).getTwistInWorldFrame().getOriginTwist();
  generalized_velocity.tail(joint_number_) = getJointVelocity();
  return generalized_velocity;
}

Eigen::VectorXd StateVariable::getGeneralizedAcceleration() const {
  const int DOF = 6;
  Eigen::VectorXd generalized_acceleration(joint_number_ + DOF);
  generalized_acceleration.head(DOF) =
      this->links_state_.at(Link::ID::kBase).getAccelInWorldFrame().getOriginAccel();
  generalized_acceleration.tail(joint_number_) = getJointAcceleration();
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

void StateVariable::setLinkTwistInWorldFrame(const int link_id, const Eigen::VectorXd &twist) {
  try {
    checkLinkIdToCall(link_id);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to set link twist.");
  }
  this->links_state_.at(link_id).setTwistInWorldFrame(Twist(Frame::kWorld, twist));
}

void StateVariable::setLinkAccelInWorldFrame(const int link_id, const Eigen::VectorXd &accel) {
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
    for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
      this->joints_state_.at(joint_id).setPosition(joint_position(joint_id));
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set joint position.");
  }
}

void StateVariable::setJointVelocity(const Eigen::VectorXd &joint_velocity) {
  try {
    checkJointVectorToSet(joint_velocity);
    for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
      this->joints_state_.at(joint_id).setVelocity(joint_velocity(joint_id));
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set joint velocity.");
  }
}

void StateVariable::setJointAcceleration(const Eigen::VectorXd &joint_acceleration) {
  try {
    checkJointVectorToSet(joint_acceleration);
    for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
      this->joints_state_.at(joint_id).setAcceleration(joint_acceleration(joint_id));
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set joint acceleration.");
  }
}

void StateVariable::setJointEffort(const Eigen::VectorXd &joint_effort) {
  try {
    checkJointVectorToSet(joint_effort);
    for (int joint_id = 0; joint_id < joint_number_; joint_id++) {
      this->joints_state_.at(joint_id).setEffort(joint_effort(joint_id));
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
    setJointVelocity(generalized_velocity.tail(joint_number_));
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
    setJointAcceleration(generalized_acceleration.tail(joint_number_));
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::invalid_argument("Error: Failed to set generalized acceleration.");
  }
}
} // namespace spacedyn_ros
