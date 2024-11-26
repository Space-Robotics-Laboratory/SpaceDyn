#include "spacedyn_ros/robot/robot.hpp"
#include "spacedyn_ros/motion/dynamics.hpp"
#include "spacedyn_ros/motion/integral.hpp"
#include "spacedyn_ros/motion/kinematics.hpp"
#include "spacedyn_ros/robot/model.hpp"
#include "spacedyn_ros/robot/state_variable.hpp"
#include <eigen3/Eigen/Core>
#include <iostream>

namespace spacedyn_ros {

Robot::Robot(Model const &model) : model_(model), state_variable_(model) {}
Robot::Robot(Model const &model, StateVariable const &state_variable)
    : model_(model), state_variable_(state_variable) {
  try {
    checkStateVariable(state_variable);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to create Robot object");
  }
}
Robot::Robot(const Robot &robot, const StateVariable &state_variable)
    : model_(robot.getModel()), state_variable_(state_variable) {
  try {
    checkStateVariable(state_variable);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to create Robot object");
  }
}

Robot::Robot(const std::string &path_to_urdf)
    : model_(Model(path_to_urdf)), state_variable_(model_) {}

void Robot::checkStateVariable(const StateVariable &state_variable) const {
  if (state_variable.getLinkNumber() != model_.getLinkage().getLinkNumber()) {
    throw std::invalid_argument("Error: Link number in StateVariable is not equal to "
                                "link number in Model");
  }
}

Eigen::MatrixXd Robot::computeGeneralizedJacobianForLink(const int link_id) const {
  return Kinematics::computeGeneralizedJacobianForLink(*this, link_id);
}

Eigen::MatrixXd Robot::computeGeneralizedJacobianForEndEffector(const int end_effector_id) const {
  return Kinematics::computeGeneralizedJacobianForEndEffector(*this, end_effector_id);
}

Eigen::MatrixXd Robot::computeGeneralizedJacobianForEndTip(const int end_effector_id) const {
  return Kinematics::computeGeneralizedJacobianForEndTip(*this, end_effector_id);
}

Eigen::MatrixXd Robot::computeJointToLinkJacobian(const int link_id) const {
  // TODO: Check if this style is suitable
  return Kinematics::computeJointToLinkJacobian(*this, link_id);
}

Eigen::MatrixXd Robot::computeJointToEndTipJacobian(const int end_effector_id) const {
  return Kinematics::computeJointToEndTipJacobian(*this, end_effector_id);
}

Eigen::MatrixXd Robot::computeBaseToLinkJacobian(const int link_id) const {
  return Kinematics::computeBaseToLinkJacobian(*this, link_id);
}

Eigen::MatrixXd Robot::computeBaseToEndTipJacobian(const int end_effector_id) const {
  return Kinematics::computeBaseToEndTipJacobian(*this, end_effector_id);
}

Eigen::MatrixXd Robot::computeInertiaMatrix() const {
  return Dynamics::computeRobotInertiaMatrix(*this);
}

Eigen::MatrixXd Robot::computeGeneralizedInertiaMatrix() const {
  return Dynamics::computeRobotGeneralizedInertiaMatrix(*this);
}

Eigen::MatrixXd Robot::computeInertiaMatrixForBaseMotion() const {
  return Dynamics::computeInertiaMatrixForBaseMotion(*this);
}

Eigen::MatrixXd Robot::computeCouplingInertiaMatrix() const {
  return Dynamics::computeCouplingInertiaMatrix(*this);
}

Eigen::MatrixXd Robot::computeInertiaMatrixForJointMotion() const {
  return Dynamics::computeInertiaMatrixForJointMotion(*this);
}

Eigen::VectorXd Robot::computeNonlinearVelocityTerm() const {
  return Dynamics::computeNonlinearVelocityTerm(*this);
}

Eigen::VectorXd Robot::computeGeneralizedNonlinearVelocityTerm() const {
  return Dynamics::computeGeneralizedNonlinearVelocityTerm(*this);
}

Eigen::VectorXd Robot::computeGravityTerm() const { return Dynamics::computeGravityTerm(*this); }

Eigen::Vector3d Robot::computeCenterOfMassInWorldFrame() const {
  return Dynamics::computeCenterOfMassInWorldFrame(*this);
}

Eigen::Vector3d Robot::computeVelocityOfCenterOfMassInWorldFrame() const {
  return Dynamics::computeVelocityOfCenterOfMassInWorldFrame(*this);
}

Eigen::Vector3d Robot::computeAccelerationOfCenterOfMassInWorldFrame() const {
  return Dynamics::computeAccelerationOfCenterOfMassInWorldFrame(*this);
}

Eigen::Vector6d Robot::computeMomentumInWorldFrame() const {
  return Dynamics::computeRobotMomentumInWorldFrame(*this);
}

Eigen::Vector6d Robot::computeMomentumAroundBaseInWorldFrame() const {
  return Dynamics::computeRobotMomentumAroundBaseInWorldFrame(*this);
}

double Robot::computeKineticEnergy() const { return Dynamics::computeRobotKineticEnergy(*this); }

StateVariable Robot::computeForward() const {
  auto sv_fk = Kinematics::computeForward(*this, true, true, false);
  auto sv_fd = Dynamics::computeForward(Robot(*this, sv_fk));
  return Kinematics::computeForward(Robot(*this, sv_fd), false, false, true);
}

Eigen::VectorXd
Robot::computeInverseDynamics(const Eigen::Vector6d &desired_base_acceleration,
                              const Eigen::VectorXd &desired_joint_acceleration) const {
  auto robot_cpy = *this;
  robot_cpy.overwriteBaseAccelInWorldFrame(desired_base_acceleration);
  robot_cpy.overwriteJointAcceleration(desired_joint_acceleration);
  robot_cpy.updateKinematics(true, true, true);
  auto sv = Dynamics::computeInverse(robot_cpy);
  auto base_ext_wrench = sv.getLinkState(Link::ID::kBase).getExternallyAppliedWrenchInWorldFrame();
  auto joint_effort = sv.getJointEffort();
  Eigen::VectorXd output(getDof());
  output.head(6) = base_ext_wrench.getWrench();
  output.tail(getActuatorNumber()) = joint_effort;
  return output;
};

Eigen::VectorXd
Robot::computeInverseDynamicsInJointSpace(const Eigen::VectorXd &desired_joint_acceleration) const {
  auto robot_cpy = *this;
  robot_cpy.overwriteJointAcceleration(desired_joint_acceleration);
  return Dynamics::computeInverseInJointSpace(robot_cpy);
}

const Model &Robot::getModel() const { return model_; }
double Robot::getTotalMass() const { return model_.getTotalMass(); }
const StateVariable &Robot::getStateVariable() const { return state_variable_; }

int Robot::getDof() const { return model_.getDof(); }
int Robot::getLinkNumber() const { return model_.getLinkNumber(); }
int Robot::getJointNumber() const { return model_.getJointNumber(); }
int Robot::getActuatorNumber() const { return model_.getActuatorNumber(); }
int Robot::getEndEffectorNumber() const { return model_.getEndEffectorNumber(); }
const Link &Robot::getLink(const int link_id) const { return model_.getLinkage().getLink(link_id); }
const Link &Robot::getEndEffector(const int end_effector_id) const {
  return model_.getLinkage().getEndEffector(end_effector_id);
}
const Joint &Robot::getJoint(const int joint_id) const {
  return model_.getLinkage().getJoint(joint_id);
}

const LinkState &Robot::getLinkState(const int link_id) const {
  return state_variable_.getLinkState(link_id);
}

const LinkState &Robot::getEndEffectorState(const int end_effector_id) const {
  return getLinkState(getEndEffector(end_effector_id).getId());
}

const JointState &Robot::getJointState(const int joint_id) const {
  return state_variable_.getJointState(joint_id);
}

Eigen::VectorXd Robot::getJointPosition() const { return state_variable_.getJointPosition(); }

Eigen::VectorXd Robot::getJointVelocity() const { return state_variable_.getJointVelocity(); }

Eigen::VectorXd Robot::getJointAcceleration() const {
  return state_variable_.getJointAcceleration();
}

Eigen::VectorXd Robot::getJointEffort() const { return state_variable_.getJointEffort(); }

Eigen::VectorXd Robot::getGeneralizedVelocity() const {
  return state_variable_.getGeneralizedVelocity();
}

Eigen::VectorXd Robot::getGeneralizedAcceleration() const {
  return state_variable_.getGeneralizedAcceleration();
}

Eigen::VectorXd Robot::getGeneralizedForce() const {
  const int CARTESIAN_DIM = 6;
  const int actuator_number = getActuatorNumber();
  Eigen::VectorXd generalized_force(getDof());
  generalized_force.head(CARTESIAN_DIM) =
      getLinkState(Link::ID::kBase).getExternallyAppliedWrenchInWorldFrame().getWrench();
  generalized_force.tail(actuator_number) = state_variable_.getJointEffort();

  for (int link_id = 1; link_id < getLinkNumber(); ++link_id) {
    auto jacob_bi = computeBaseToLinkJacobian(link_id);
    auto jacob_mi = computeJointToLinkJacobian(link_id);
    Eigen::MatrixXd jacobian_trans(getDof(), CARTESIAN_DIM);
    jacobian_trans.topRows(CARTESIAN_DIM) = jacob_bi.transpose();
    jacobian_trans.bottomRows(actuator_number) = jacob_mi.transpose();
    Eigen::VectorXd link_ext_force =
        getLinkState(link_id).getExternallyAppliedWrenchInWorldFrame().getWrench();
    generalized_force += jacobian_trans * link_ext_force;
  }

  return generalized_force;
}

const Eigen::Isometry3d &Robot::getBasePose() const {
  return getLinkState(Link::ID::kBase).getPoseInWorldFrame().getPoseInWorldFrame();
}

const Pose &Robot::getBasePoseInWorldFrame() const {
  return getLinkState(Link::ID::kBase).getPoseInWorldFrame();
}

const Pose &Robot::getLinkPoseInWorldFrame(const int link_id) const {
  return getLinkState(link_id).getPoseInWorldFrame();
}

const Pose &Robot::getEndEffectorPoseInWorldFrame(const int end_effector_id) const {
  return getEndEffectorState(end_effector_id).getPoseInWorldFrame();
}

const Pose Robot::getEndTipPoseInWorldFrame(const int end_effector_id) const {
  auto ee_pose = getEndEffectorPoseInWorldFrame(end_effector_id);
  auto tf_to_end_tip = getEndEffector(end_effector_id).getTransformToEndTip();
  return ee_pose.computePointPose(tf_to_end_tip);
}

const Twist &Robot::getBaseTwistInWorldFrame() const {
  return getLinkState(Link::ID::kBase).getTwistInWorldFrame();
}

Twist Robot::getBaseTwistInLocalFrame() const {
  return getLinkState(Link::ID::kBase).getTwistInLocalFrame();
}

const Twist &Robot::getLinkTwistInWorldFrame(const int link_id) const {
  return getLinkState(link_id).getTwistInWorldFrame();
}

Twist Robot::getLinkTwistInLocalFrame(const int link_id) const {
  return getLinkState(link_id).getTwistInLocalFrame();
}

const Twist &Robot::getEndEffectorTwistInWorldFrame(const int end_effector_id) const {
  return getEndEffectorState(end_effector_id).getTwistInWorldFrame();
}

Twist Robot::getEndEffectorTwistInLocalFrame(const int end_effector_id) const {
  return getEndEffectorState(end_effector_id).getTwistInLocalFrame();
}

Twist Robot::getEndTipTwistInWorldFrame(const int end_effector_id) const {
  // TODO: make it more efficient
  auto ee_twist = getEndEffectorTwistInWorldFrame(end_effector_id);
  auto ee_pose = getEndEffectorPoseInWorldFrame(end_effector_id);
  auto tip_pos = getEndTipPoseInWorldFrame(end_effector_id);
  return ee_twist.computePointTwist(
      ee_pose.computeTransformToPoint(Frame::kWorld, tip_pos).getTranslation());
}

const Accel &Robot::getBaseAccelInWorldFrame() const {
  return getLinkState(Link::ID::kBase).getAccelInWorldFrame();
}

Accel Robot::getBaseAccelInLocalFrame() const {
  return getLinkState(Link::ID::kBase).getAccelInLocalFrame();
}

const Accel &Robot::getLinkAccelInWorldFrame(const int link_id) const {
  return getLinkState(link_id).getAccelInWorldFrame();
}

Accel Robot::getLinkAccelInLocalFrame(const int link_id) const {
  return getLinkState(link_id).getAccelInLocalFrame();
}

const Accel &Robot::getEndEffectorAccelInWorldFrame(const int end_effector_id) const {
  return getEndEffectorState(end_effector_id).getAccelInWorldFrame();
}

Accel Robot::getEndEffectorAccelInLocalFrame(const int end_effector_id) const {
  return getEndEffectorState(end_effector_id).getAccelInLocalFrame();
}

Accel Robot::getEndTipAccelInWorldFrame(const int end_effector_id) const {
  auto ee_accel = getEndEffectorAccelInWorldFrame(end_effector_id);
  auto ee_twist = getEndEffectorTwistInWorldFrame(end_effector_id);
  auto tf_to_end_tip = getEndEffector(end_effector_id).getTransformToEndTip();
  return ee_accel.computePointAccel(ee_twist, tf_to_end_tip.getTransform().translation());
}

void Robot::step() {
  // TODO: Add system to choose integration method
  auto sv = Integral::rungeKutta4(*this);
  // auto sv = Integral::euler(*this);

  setStateVariable(sv);
  clearAllLinkExternallyAppliedWrench();
  clearJointEffort();
}

void Robot::updateKinematics(bool pose, bool twist, bool accel) {
  state_variable_ = Kinematics::computeForward(*this, pose, twist, accel);
}

void Robot::setStateVariable(const StateVariable &state_variable) {
  try {
    checkStateVariable(state_variable);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    throw std::runtime_error("Error: Failed to set state variable");
  }

  this->state_variable_ = state_variable;
}

void Robot::overwriteBaseState(const LinkState &base_state) {
  state_variable_.setLinkState(Link::ID::kBase, base_state);
}

void Robot::overwriteBasePoseInWorldFrame(const Eigen::Isometry3d &base_pose) {
  state_variable_.setLinkPoseInWorldFrame(Link::ID::kBase, base_pose);
}

void Robot::overwriteBasePoseInWorldFrame(const Eigen::Vector3d &base_position,
                                          const Eigen::Matrix3d &base_attitude) {
  state_variable_.setLinkPoseInWorldFrame(
      Link::ID::kBase,
      Eigen::Isometry3d(Eigen::Translation3d(base_position) * Eigen::Quaterniond(base_attitude)));
}

void Robot::overwriteBaseTwistInWorldFrame(const Eigen::Vector6d &base_twist) {
  state_variable_.setLinkTwistInWorldFrame(Link::ID::kBase, base_twist);
}

void Robot::overwriteBaseAccelInWorldFrame(const Eigen::Vector6d &base_accel) {
  state_variable_.setLinkAccelInWorldFrame(Link::ID::kBase, base_accel);
}

void Robot::overwriteLinkPoseInWorldFrame(const int link_id, const Eigen::Isometry3d &link_pose) {
  state_variable_.setLinkPoseInWorldFrame(link_id, link_pose);
}

void Robot::overwriteLinkTwistInWorldFrame(const int link_id, const Eigen::Vector6d &link_twist) {
  state_variable_.setLinkTwistInWorldFrame(link_id, link_twist);
}

void Robot::overwriteLinkAccelInWorldFrame(const int link_id, const Eigen::Vector6d &link_accel) {
  state_variable_.setLinkAccelInWorldFrame(link_id, link_accel);
}

void Robot::overwriteJointPosition(const Eigen::VectorXd &joint_position) {
  state_variable_.setJointPosition(joint_position);
}

void Robot::overwriteJointVelocity(const Eigen::VectorXd &joint_velocity) {
  state_variable_.setJointVelocity(joint_velocity);
}

void Robot::overwriteJointAcceleration(const Eigen::VectorXd &joint_acceleration) {
  state_variable_.setJointAcceleration(joint_acceleration);
}

void Robot::clearBaseState() { clearLinkState(Link::ID::kBase); }

void Robot::clearBaseTwist() { clearLinkTwist(Link::ID::kBase); }

void Robot::clearBaseAccel() { clearLinkAccel(Link::ID::kBase); }

void Robot::clearBaseExternallyAppliedWrench() {
  clearLinkExternallyAppliedWrench(Link::ID::kBase);
}

void Robot::clearLinkState(const int link_id) {
  state_variable_.setLinkState(link_id, LinkState());
}

void Robot::clearLinkTwist(const int link_id) {
  state_variable_.setLinkTwistInWorldFrame(link_id, Eigen::Vector6d::Zero());
}

void Robot::clearLinkAccel(const int link_id) {
  state_variable_.setLinkAccelInWorldFrame(link_id, Eigen::Vector6d::Zero());
}

void Robot::clearLinkExternallyAppliedWrench(const int link_id) {
  state_variable_.setLinkExternallyAppliedWrenchInWorldFrame(link_id, Eigen::Vector6d::Zero());
}

void Robot::clearAllLinkAccel() {
  for (int link_id = 0; link_id < getLinkNumber(); ++link_id) {
    clearLinkAccel(link_id);
  }
}

void Robot::clearAllLinkExternallyAppliedWrench() {
  for (int link_id = 0; link_id < getLinkNumber(); ++link_id) {
    clearLinkExternallyAppliedWrench(link_id);
  }
}

void Robot::clearJointPosition() {
  state_variable_.setJointPosition(Eigen::VectorXd::Zero(getActuatorNumber()));
}

void Robot::clearJointVelocity() {
  state_variable_.setJointVelocity(Eigen::VectorXd::Zero(getActuatorNumber()));
}

void Robot::clearJointAcceleration() {
  state_variable_.setJointAcceleration(Eigen::VectorXd::Zero(getActuatorNumber()));
}

void Robot::clearJointEffort() {
  state_variable_.setJointEffort(Eigen::VectorXd::Zero(getActuatorNumber()));
}

void Robot::applyExternalWrench(const int link_id, const Wrench &external_wrench) {
  LinkState link_state = getLinkState(link_id);
  link_state.setExternallyAppliedWrenchInWorldFrame(external_wrench);
  state_variable_.setLinkState(link_id, link_state);
}

void Robot::applyJointEffort(const Eigen::VectorXd &joint_effort) {
  state_variable_.setJointEffort(joint_effort);
}

geometry_msgs::msg::TransformStamped Robot::basePoseToRosTf() const {
  geometry_msgs::msg::TransformStamped msg;
  auto base = getLink(Link::ID::kBase);
  auto base_pose = getBasePoseInWorldFrame();
  msg = base_pose.toRosMessage("world", base.getName());
  return msg;
}

std::vector<geometry_msgs::msg::TransformStamped> Robot::jointPoseToRosTf() const {
  // TODO: Add local transform. Currently, all joints are in world frame.
  const int joint_number = getJointNumber();
  std::vector<geometry_msgs::msg::TransformStamped> msg(joint_number);
  auto parent_name = getLink(Link::ID::kBase).getName();
  for (int joint_id = 0; joint_id < joint_number; ++joint_id) {
    auto joint = getJoint(joint_id);
    auto joint_pose = getJointState(joint_id).getPoseInWorldFrame();
    msg[joint_id] = joint_pose.toRosMessage("world", joint.getName());
  }
  return msg;
}

void Robot::overWriteJointState(const sensor_msgs::msg::JointState &joint_state_msg) {
  // TODO: Check id joint_state is valid. Move this to children.
  Eigen::VectorXd joint_position(joint_state_msg.position.size());
  Eigen::VectorXd joint_velocity(joint_state_msg.velocity.size());
  Eigen::VectorXd joint_acceleration(joint_state_msg.effort.size());
  for (int i = 0; i < static_cast<int>(joint_state_msg.position.size()); i++) {
    joint_position(i) = joint_state_msg.position[i];
    joint_velocity(i) = joint_state_msg.velocity[i];
    joint_acceleration(i) = joint_state_msg.effort[i];
  }
  overwriteJointPosition(joint_position);
  overwriteJointVelocity(joint_velocity);
  overwriteJointAcceleration(joint_acceleration);
}

} // namespace spacedyn_ros
