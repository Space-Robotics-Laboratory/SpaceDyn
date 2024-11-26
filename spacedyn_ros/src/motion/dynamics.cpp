#include "spacedyn_ros/motion/dynamics.hpp"
#include "spacedyn_ros/motion/kinematics.hpp"
#include "spacedyn_ros/util/matrix_operation.hpp"
#include <iostream>

namespace spacedyn_ros {
Dynamics::Dynamics() {}

// TODO: For the sake of speed, consider (m, sv) instead of (robot)
LinkState Dynamics::inverseLinkWrench(const Link &link, const LinkState &link_state) {
  Inertia link_inertia = link.computeInertiaInWorldFrame(link_state.getPoseInWorldFrame());
  Accel link_accel = link_state.getAccelInWorldFrame();
  Twist link_twist = link_state.getTwistInWorldFrame();
  Eigen::Vector6d wrench;
  // Compute wrench from Newton-Euler equation
  wrench.head(3) = link_inertia.getMass() * (link_accel.getLinearAcceleration()); // F = ma
  wrench.tail(3) =
      link_inertia.getInertiaTensor() * link_accel.getAngularAcceleration() +
      link_twist.getAngularVelocity().cross(link_inertia.getInertiaTensor() *
                                            link_twist.getAngularVelocity()); // T = I*dw + w x I*w
  LinkState output = link_state;
  output.setTotalWrenchOnLinkInWorldFrame(Wrench(Frame::kWorld, wrench));
  return output;
}

JointState Dynamics::inverseJointWrench(const Joint &joint, const JointState &joint_state,
                                        const Link &link, const LinkState &link_state,
                                        const Wrench &summed_wrench_by_children_joints_at_link,
                                        const Eigen::Vector3d &gravity) {
  Eigen::Vector3d trans_from_parent_joint_to_link_com =
      joint_state.getPoseInWorldFrame().computeTranslationToPoint(
          link_state.getPoseInWorldFrame().getPoseInWorldFrame());

  // Compute the required wrench at CoM applied by the joint.
  // Wrench at link CoM: (link total) = (joint)      + (children joints) + (gravity) + (external)
  //                 <=> (joint)      = (link total) - (children joints) - (gravity) - (external)
  Wrench gravity_wrench(Frame::kWorld, link.getInertiaInLocalFrame().getMass() * gravity,
                        Eigen::Vector3d::Zero());
  // TODO: Check Joint wrench and effort direction
  Wrench joint_wrench_at_link_com = link_state.getTotalWrenchOnLinkInWorldFrame() -
                                    summed_wrench_by_children_joints_at_link - gravity_wrench -
                                    link_state.getExternallyAppliedWrenchInWorldFrame();

  // Ta = Pab x Fb
  Wrench joint_wrench = joint_wrench_at_link_com.computeWrenchByInvertingPointWrench(
      trans_from_parent_joint_to_link_com);

  JointState output = joint_state;
  output.setWrenchToChildInWorldFrame(joint_wrench);
  double effort = joint.computeEffortToAchieveWrench(output);
  output.setEffort(effort);
  return output;
}

LinkState Dynamics::inverseBaseExtWrench(const Link &link, const LinkState &link_state,
                                         const Wrench &summed_wrench_by_children_joints_at_link,
                                         const Eigen::Vector3d &gravity) {
  if (link.getId() != Link::ID::kBase) {
    throw std::invalid_argument("Error: Link is not a base link.");
  }
  Wrench gravity_wrench(Frame::kWorld, link.getInertiaInLocalFrame().getMass() * gravity,
                        Eigen::Vector3d::Zero());

  Wrench ext_wrench_at_base_com = link_state.getTotalWrenchOnLinkInWorldFrame() -
                                  summed_wrench_by_children_joints_at_link - gravity_wrench;

  LinkState output = link_state;
  output.setExternallyAppliedWrenchInWorldFrame(ext_wrench_at_base_com);
  return output;
}

StateVariable Dynamics::computeForward(const Robot &robot) {
  // Compute forward using Equation of Motion derived from Lagrangian
  // H * q_ddot + C(q, q_dot) + g(q) = F
  // F = (joint effort) + (sum of external wrench) => Generalized force
  StateVariable result_sv = robot.getStateVariable();
  auto H = computeRobotInertiaMatrix(robot);

  // For non-linear and gravity term, use inverse dynamics
  // (recursive Newton-Euler): C(q, q_dot) + g(q) = F' (q_ddot = 0, F_ext = 0)
  // TODO: Rewrite this part using Robot.computeInverse()
  Robot robot_no_acc = robot;
  robot_no_acc.clearBaseAccel();         // q_ddot = 0
  robot_no_acc.clearJointAcceleration(); // q_ddot = 0
  robot_no_acc.setStateVariable(Kinematics::computeForward(robot_no_acc, false, false, true));
  robot_no_acc.clearAllLinkExternallyAppliedWrench(); // F_ext = 0
  robot_no_acc.setStateVariable(computeInverse(robot_no_acc));
  Eigen::VectorXd non_linear_plus_gravity_term = robot_no_acc.getGeneralizedForce();

  Eigen::VectorXd inertia_force = robot.getGeneralizedForce() - non_linear_plus_gravity_term;
  Eigen::VectorXd generalized_acceleration = H.inverse() * inertia_force;

  // TODO: Replace this part with setGeneralizedAcceleration()
  const int SPACE_DIM = 6;
  result_sv.setLinkAccelInWorldFrame(Link::ID::kBase, generalized_acceleration.head(SPACE_DIM));
  result_sv.setJointAcceleration(generalized_acceleration.tail(robot.getActuatorNumber()));

  return result_sv;
}

StateVariable Dynamics::computeInverse(const Robot &robot) {
  StateVariable result_sv = robot.getStateVariable();
  Eigen::Vector3d gravity = robot.getModel().getGravity();

  // Start recursive computation from the first end-effector
  int end_effector_id = 0;
  Link link = robot.getEndEffector(end_effector_id);
  int loop_count = 0;
  std::vector<Wrench> summed_wrench_by_children_joints_at_link(robot.getLinkNumber());

  // stock children number to decide jump or not
  std::vector<int> link_jump_count(robot.getLinkNumber());
  for (int i = 0; i < robot.getLinkNumber(); i++) {
    link_jump_count.at(i) = robot.getLink(i).getChildrenNumber();
  }

  // Skip base computation because no joint on base
  // TODO: Make it better by using Connection class
  while (link.getId() != Link::ID::kBase) {
    int link_id = link.getId();
    int joint_id = link.getParentJointId();
    Joint joint = robot.getJoint(joint_id);
    LinkState link_state = result_sv.getLinkState(link_id);
    JointState joint_state = result_sv.getJointState(joint_id);
    LinkState link_state_idyn_done = inverseLinkWrench(link, link_state);

    // Compute inverse dynamics on link and joint
    JointState joint_state_idyn_done =
        inverseJointWrench(joint, joint_state, link, link_state_idyn_done,
                           summed_wrench_by_children_joints_at_link.at(link_id), gravity);

    // Set the output state
    result_sv.setLinkState(link_id, link_state_idyn_done);
    result_sv.setJointState(joint_id, joint_state_idyn_done);

    // Handle link order to compute and sum-up children joints wrench
    int parent_link_id = link.getParentId();
    Link parent_link = robot.getLink(parent_link_id);
    LinkState parent_link_state = result_sv.getLinkState(parent_link_id);

    Eigen::Vector3d trans_parent_link_to_joint =
        parent_link_state.getPoseInWorldFrame().computeTranslationToPoint(
            joint_state_idyn_done.getPoseInWorldFrame());

    // Add inverse of joint wrench to the parent link: Newton 3rd law
    // Sum of them is wrench from all children
    auto wrench_by_child_joint = -joint_state_idyn_done.getWrenchToChildInWorldFrame();
    summed_wrench_by_children_joints_at_link.at(parent_link_id) +=
        wrench_by_child_joint.computeWrenchByInvertingPointWrench(trans_parent_link_to_joint);

    // If parent link has other child, stash the link and jump to the next end effector to compute
    // summed wrench by children before compute the parent wrench.
    // So, all children should have been computed to computed parent.
    bool parent_has_uncomputed_child = link_jump_count.at(parent_link_id) > 1;
    if (parent_has_uncomputed_child) {
      // Stash the link and jump to the next end effector
      link_jump_count.at(parent_link_id)--;
      end_effector_id++;
      link = robot.getEndEffector(end_effector_id);
    } else {
      // Simply got parent link without jumping
      link = parent_link;
    }

    // For safety to avoid infinite loop
    loop_count++;
    if (loop_count >= result_sv.getLinkNumber()) {
      throw std::logic_error("Error: Unexpected error. Loop count is over link number");
      break;
    }
  }

  // Compute base link's external wrench
  // link == base here
  LinkState link_state_link_idyn_done = result_sv.getLinkState(Link::ID::kBase);
  LinkState link_state_base_idyn_done = inverseLinkWrench(link, link_state_link_idyn_done);
  LinkState link_state_base_external_idyn_done =
      inverseBaseExtWrench(link, link_state_base_idyn_done,
                           summed_wrench_by_children_joints_at_link.at(Link::ID::kBase), gravity);
  result_sv.setLinkState(Link::ID::kBase, link_state_base_external_idyn_done);
  return result_sv;
}

Eigen::VectorXd Dynamics::computeInverseInJointSpace(const Robot &robot) {
  // TODO: External wrench is not considered in this function currently
  auto GH = computeRobotGeneralizedInertiaMatrix(robot);
  auto GC = computeGeneralizedNonlinearVelocityTerm(robot);
  auto GG = computeGeneralizedGravityTerm(robot);
  auto des_acc = robot.getJointAcceleration();
  auto des_effort = GH * des_acc + GC + GG;
  return des_effort;
}

Eigen::Vector3d Dynamics::computeCenterOfMassInWorldFrame(const Robot &robot) {
  Eigen::Vector3d com = Eigen::Vector3d::Zero();
  double total_mass = robot.getTotalMass();
  for (int i = 0; i < robot.getLinkNumber(); i++) {
    Link link = robot.getLink(i);
    LinkState link_state = robot.getLinkState(i);
    double mass = link.getInertiaInLocalFrame().getMass();
    com += link_state.getPoseInWorldFrame().getPositionInWorldFrame() * mass;
  }
  return com / total_mass;
}

Eigen::Vector3d Dynamics::computeVelocityOfCenterOfMassInWorldFrame(const Robot &robot) {
  Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
  double total_mass = robot.getTotalMass();
  for (int i = 0; i < robot.getLinkNumber(); i++) {
    Link link = robot.getLink(i);
    LinkState link_state = robot.getLinkState(i);
    double mass = link.getInertiaInLocalFrame().getMass();
    velocity += link_state.getTwistInWorldFrame().getLinearVelocity() * mass;
  }
  return velocity / total_mass;
}

Eigen::Vector3d Dynamics::computeAccelerationOfCenterOfMassInWorldFrame(const Robot &robot) {
  Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();
  double total_mass = robot.getTotalMass();
  for (int i = 0; i < robot.getLinkNumber(); i++) {
    Link link = robot.getLink(i);
    LinkState link_state = robot.getLinkState(i);
    double mass = link.getInertiaInLocalFrame().getMass();
    acceleration += link_state.getAccelInWorldFrame().getLinearAcceleration() * mass;
  }
  return acceleration / total_mass;
}

Eigen::MatrixXd Dynamics::computeInertiaMatrixForBaseMotion(const Robot &robot) {
  Eigen::MatrixXd inertia_matrix = Eigen::MatrixXd::Zero(6, 6);
  double total_mass = robot.getTotalMass();
  Pose base_pose = robot.getBasePoseInWorldFrame();
  Eigen::Matrix3d base_to_com_skew_symmetric_matrix =
      skewSymmetric(computeCenterOfMassInWorldFrame(robot) - base_pose.getPositionInWorldFrame());

  Eigen::Matrix3d link_inertia_to_base_rot = Eigen::Matrix3d::Zero();

  for (int link_id = 0; link_id < robot.getLinkNumber(); link_id++) {
    Link link = robot.getLink(link_id);
    Pose link_pose = robot.getLinkState(link_id).getPoseInWorldFrame();
    Eigen::Matrix3d base_to_link_skew_symmetric_matrix =
        skewSymmetric(link_pose.getPositionInWorldFrame() - base_pose.getPositionInWorldFrame());

    link_inertia_to_base_rot +=
        link.computeInertiaInWorldFrame(link_pose).getInertiaTensor() -
        link.getMass() * base_to_link_skew_symmetric_matrix * base_to_link_skew_symmetric_matrix;
  }

  inertia_matrix.topLeftCorner(3, 3) = total_mass * Eigen::Matrix3d::Identity();
  inertia_matrix.bottomRightCorner(3, 3) = link_inertia_to_base_rot;
  inertia_matrix.topRightCorner(3, 3) = -total_mass * base_to_com_skew_symmetric_matrix;
  inertia_matrix.bottomLeftCorner(3, 3) = total_mass * base_to_com_skew_symmetric_matrix;

  return inertia_matrix;
}

Eigen::MatrixXd Dynamics::computeInertiaMatrixForJointMotion(const Robot &robot) {
  Eigen::MatrixXd inertia_matrix =
      Eigen::MatrixXd::Zero(robot.getActuatorNumber(), robot.getActuatorNumber());
  for (int link_id = 0; link_id < robot.getLinkNumber(); link_id++) {
    Link link = robot.getLink(link_id);
    Eigen::MatrixXd jacobian = robot.computeJointToLinkJacobian(link_id);
    inertia_matrix += jacobian.bottomRows(3).transpose() *
                          link.computeInertiaInWorldFrame(robot.getLinkPoseInWorldFrame(link_id))
                              .getInertiaTensor() *
                          jacobian.bottomRows(3) +
                      link.getMass() * jacobian.topRows(3).transpose() * jacobian.topRows(3);
  }
  return inertia_matrix;
}

Eigen::MatrixXd Dynamics::computeCouplingInertiaMatrix(const Robot &robot) {
  const int CARTESIAN_DIM = 6;
  Eigen::MatrixXd coupling_inertia_matrix =
      Eigen::MatrixXd::Zero(CARTESIAN_DIM, robot.getActuatorNumber());
  Pose base_pose = robot.getBasePoseInWorldFrame();

  for (int link_id = 0; link_id < robot.getLinkNumber(); link_id++) {
    Eigen::MatrixXd coupling_inertia_matrix_link =
        Eigen::MatrixXd::Zero(CARTESIAN_DIM, robot.getActuatorNumber());

    Link link = robot.getLink(link_id);
    Pose link_pose = robot.getLinkState(link_id).getPoseInWorldFrame();
    Eigen::MatrixXd jacobian = robot.computeJointToLinkJacobian(link_id);

    Eigen::Matrix3d base_to_link_skew_symmetric_matrix =
        skewSymmetric(link_pose.getPositionInWorldFrame() - base_pose.getPositionInWorldFrame());

    coupling_inertia_matrix_link.topRows(3) = link.getMass() * jacobian.topRows(3);
    coupling_inertia_matrix_link.bottomRows(3) =
        link.computeInertiaInWorldFrame(link_pose).getInertiaTensor() * jacobian.bottomRows(3) +
        link.getMass() * base_to_link_skew_symmetric_matrix * jacobian.topRows(3);

    coupling_inertia_matrix += coupling_inertia_matrix_link;
  }
  return coupling_inertia_matrix;
}

Eigen::MatrixXd Dynamics::computeRobotInertiaMatrix(const Robot &robot) {
  const int CARTESIAN_DIM = 6;
  const int actuator_number = robot.getActuatorNumber();
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(robot.getDof(), robot.getDof());
  auto Hb = computeInertiaMatrixForBaseMotion(robot);
  auto Hm = computeInertiaMatrixForJointMotion(robot);
  auto Hbm = computeCouplingInertiaMatrix(robot);

  H.topLeftCorner(CARTESIAN_DIM, CARTESIAN_DIM) = Hb;
  H.bottomRightCorner(actuator_number, actuator_number) = Hm;
  H.topRightCorner(CARTESIAN_DIM, actuator_number) = Hbm;
  H.bottomLeftCorner(actuator_number, CARTESIAN_DIM) = Hbm.transpose();
  return H;
}

Eigen::MatrixXd Dynamics::computeRobotGeneralizedInertiaMatrix(const Robot &robot) {
  auto Hb = computeInertiaMatrixForBaseMotion(robot);
  auto Hm = computeInertiaMatrixForJointMotion(robot);
  auto Hbm = computeCouplingInertiaMatrix(robot);
  return Hm - Hbm.transpose() * Hb.inverse() * Hbm;
}

Eigen::VectorXd Dynamics::computeNonlinearVelocityTerm(const Robot &robot) {
  // Compute the non-linear velocity term using the inverse dynamics
  // C(q, q_dot) = F' (q_ddot = 0, F_ext = 0, g = 0)
  Model zeroG_model = robot.getModel();
  zeroG_model.setGravity(Eigen::Vector3d::Zero()); // g = 0
  Robot robot_cpy(zeroG_model);
  robot_cpy.setStateVariable(robot.getStateVariable());
  robot_cpy.clearBaseAccel();                      // q_ddot = 0
  robot_cpy.clearJointAcceleration();              // q_ddot = 0
  robot_cpy.clearAllLinkExternallyAppliedWrench(); // F_ext = 0
  robot_cpy.updateKinematics(true, true, true);
  robot_cpy.setStateVariable(computeInverse(robot_cpy));
  Eigen::VectorXd generalized_force = robot_cpy.getGeneralizedForce();
  return generalized_force;
}

Eigen::VectorXd Dynamics::computeGeneralizedNonlinearVelocityTerm(const Robot &robot) {
  auto C = robot.computeNonlinearVelocityTerm();
  Eigen::Vector6d Cb = C.head(6);
  Eigen::VectorXd Cm = C.tail(robot.getActuatorNumber());
  auto Hb = computeInertiaMatrixForBaseMotion(robot);
  auto Hbm = computeCouplingInertiaMatrix(robot);

  return Cm - Hbm.transpose() * Hb.inverse() * Cb;
}

Eigen::VectorXd Dynamics::computeGravityTerm(const Robot &robot) {
  // Compute the gravity term using the inverse dynamics
  // g(q) = F' (q_ddot = 0, F_ext = 0)
  Robot robot_cpy = robot;
  robot_cpy.clearBaseTwist();                      // q_dot = 0
  robot_cpy.clearJointVelocity();                  // q_dot = 0
  robot_cpy.clearBaseAccel();                      // q_ddot = 0
  robot_cpy.clearJointAcceleration();              // q_ddot = 0
  robot_cpy.clearAllLinkExternallyAppliedWrench(); // F_ext = 0
  robot_cpy.updateKinematics(true, true, true);
  robot_cpy.setStateVariable(computeInverse(robot_cpy));
  Eigen::VectorXd generalized_force = robot_cpy.getGeneralizedForce();
  return generalized_force;
}

Eigen::VectorXd Dynamics::computeGeneralizedGravityTerm(const Robot &robot) {
  auto G = robot.computeGravityTerm();
  Eigen::Vector6d Gb = G.head(6);
  Eigen::VectorXd Gm = G.tail(robot.getActuatorNumber());
  auto Hb = computeInertiaMatrixForBaseMotion(robot);
  auto Hbm = computeCouplingInertiaMatrix(robot);

  return Gm - Hbm.transpose() * Hb.inverse() * Gb;
}

Eigen::Vector6d Dynamics::computeRobotMomentumInWorldFrame(const Robot &robot) {
  Eigen::Vector6d momentum = Eigen::Vector6d::Zero();
  for (int link_id = 0; link_id < robot.getLinkNumber(); link_id++) {
    Link link = robot.getLink(link_id);
    LinkState link_state = robot.getLinkState(link_id);
    Inertia inertia = link.computeInertiaInWorldFrame(link_state.getPoseInWorldFrame());
    Twist twist = link_state.getTwistInWorldFrame();
    Pose pose = link_state.getPoseInWorldFrame();
    Eigen::Vector3d P = inertia.getMass() * twist.getLinearVelocity();
    Eigen::Vector3d L = inertia.getInertiaTensor() * twist.getAngularVelocity();
    momentum.head(3) += P;
    momentum.tail(3) += L + pose.getPositionInWorldFrame().cross(P);
  }
  return momentum;
}

Eigen::Vector6d Dynamics::computeRobotMomentumAroundBaseInWorldFrame(const Robot &robot) {
  Eigen::Vector6d momentum = Eigen::Vector6d::Zero();
  Eigen::Vector3d r0 = robot.getBasePoseInWorldFrame().getPositionInWorldFrame();
  for (int link_id = 0; link_id < robot.getLinkNumber(); link_id++) {
    Link link = robot.getLink(link_id);
    LinkState link_state = robot.getLinkState(link_id);
    Inertia inertia = link.computeInertiaInWorldFrame(link_state.getPoseInWorldFrame());
    Twist twist = link_state.getTwistInWorldFrame();
    Eigen::Vector3d r0i = link_state.getPoseInWorldFrame().getPositionInWorldFrame() - r0;
    Eigen::Vector3d P = inertia.getMass() * twist.getLinearVelocity();
    Eigen::Vector3d L = inertia.getInertiaTensor() * twist.getAngularVelocity();
    momentum.head(3) += P;
    momentum.tail(3) += L + r0i.cross(P);
  }
  return momentum;
}

double Dynamics::computeRobotKineticEnergy(const Robot &robot) {
  StateVariable state_variable = robot.getStateVariable();
  double kinetic_energy = 0;

  // Compute kinetic energy for each link
  for (int link_id = 0; link_id < robot.getLinkNumber(); link_id++) {
    Link link = robot.getLink(link_id);
    LinkState link_state = state_variable.getLinkState(link_id);
    Inertia inertia = link.computeInertiaInWorldFrame(link_state.getPoseInWorldFrame());
    Twist twist = link_state.getTwistInWorldFrame();
    kinetic_energy += 0.5 * inertia.getMass() * twist.getLinearVelocity().squaredNorm() +
                      0.5 * twist.getAngularVelocity().dot(inertia.getInertiaTensor() *
                                                           twist.getAngularVelocity());
  }
  return kinetic_energy;
}
} // namespace spacedyn_ros
