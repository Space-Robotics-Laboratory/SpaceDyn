#include "spacedyn_ros/motion/kinematics.hpp"
#include "spacedyn_ros/linkage/link.hpp"
#include "spacedyn_ros/robot/model.hpp"
#include "spacedyn_ros/robot/robot.hpp"
#include "spacedyn_ros/robot/state_variable.hpp"
#include "spacedyn_ros/util/matrix_operation.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <iostream>

namespace spacedyn_ros {

Kinematics::Kinematics() {}

JointState Kinematics::forwardJointPose(const Joint &joint, const LinkState &parent_link_state,
                                        const JointState &joint_state) {
  Pose joint_pose = parent_link_state.getPoseInWorldFrame()
                        .computePointPose(joint.getTransformFromParentLinkCom())
                        .computePointPose(joint.transformByActuation(joint_state));

  JointState output = joint_state;
  output.setPoseInWorldFrame(joint_pose);
  return output;
}

JointState Kinematics::forwardJointTwist(const Joint &joint, const LinkState &parent_link_state,
                                         const JointState &joint_state) {
  // TODO: Implement the local frame and stop using Eigen vector directory
  Eigen::Vector3d trans_link_to_joint =
      parent_link_state.getPoseInWorldFrame().computeTranslationToPoint(
          joint_state.getPoseInWorldFrame());
  Twist joint_twist =
      parent_link_state.getTwistInWorldFrame().computePointTwist(trans_link_to_joint) +
      joint.twistByActuation(joint_state);
  JointState output = joint_state;
  output.setTwistInWorldFrame(joint_twist);
  return output;
}

JointState Kinematics::forwardJointAccel(const Joint &joint, const LinkState &parent_link_state,
                                         const JointState &joint_state) {
  Twist parent_twist = parent_link_state.getTwistInWorldFrame();
  Accel parent_accel = parent_link_state.getAccelInWorldFrame();
  Eigen::Vector3d trans_link_to_joint =
      parent_link_state.getPoseInWorldFrame().computeTranslationToPoint(
          joint_state.getPoseInWorldFrame());
  Accel joint_accel = parent_accel.computePointAccel(parent_twist, trans_link_to_joint) +
                      joint.accelByActuation(joint_state);

  JointState output = joint_state;
  output.setAccelInWorldFrame(joint_accel);
  return output;
}

LinkState Kinematics::forwardLinkPose(const Link &link, const JointState &parent_joint_state,
                                      const LinkState &link_state) {
  Pose link_pose_at_com =
      parent_joint_state.getPoseInWorldFrame().computePointPose(link.getTransformFromParentJoint());
  LinkState output = link_state;
  output.setPoseInWorldFrame(link_pose_at_com);
  return output;
}

LinkState Kinematics::forwardLinkTwist(const Link &link, const JointState &parent_joint_state,
                                       const LinkState &link_state) {
  Eigen::Vector3d trans_joint_to_link =
      parent_joint_state.getPoseInWorldFrame().computeTranslationToPoint(
          link_state.getPoseInWorldFrame());
  Twist link_twist =
      parent_joint_state.getTwistInWorldFrame().computePointTwist(trans_joint_to_link);

  LinkState output = link_state;
  output.setTwistInWorldFrame(link_twist);
  return output;
}

LinkState Kinematics::forwardLinkAccel(const Link &link, const JointState &parent_joint_state,
                                       const LinkState &link_state) {
  Twist parent_twist = parent_joint_state.getTwistInWorldFrame();
  Accel parent_accel = parent_joint_state.getAccelInWorldFrame();
  Eigen::Vector3d trans_joint_to_link =
      parent_joint_state.getPoseInWorldFrame().computeTranslationToPoint(
          link_state.getPoseInWorldFrame());
  Accel link_accel = parent_accel.computePointAccel(parent_twist, trans_joint_to_link);

  LinkState output = link_state;
  output.setAccelInWorldFrame(link_accel);
  return output;
}

// TODO: Modify a flag to skip twist and accel computation
StateVariable Kinematics::computeForward(const Robot &robot, const bool compute_pose,
                                         const bool compute_twist, const bool compute_accel) {
  // Set current state
  StateVariable state_variable_new = robot.getStateVariable();
  Linkage linkage = robot.getModel().getLinkage();

  // Compute forward kinematics
  // Do not update base (link_id = 0) state here
  for (int link_id = 1; link_id < linkage.getLinkNumber(); link_id++) {
    Link link = linkage.getLink(link_id);
    int joint_id = link.getParentJointId();
    Joint joint = linkage.getJoint(joint_id);
    LinkState link_state = state_variable_new.getLinkState(link_id);
    JointState joint_state = state_variable_new.getJointState(joint_id);

    LinkState parent_link_state_fk_done = state_variable_new.getLinkState(link.getParentId());

    // Compute pose
    JointState joint_state_pos_done;
    LinkState link_state_pos_done;
    if (compute_pose) {
      joint_state_pos_done = forwardJointPose(joint, parent_link_state_fk_done, joint_state);
      link_state_pos_done = forwardLinkPose(link, joint_state_pos_done, link_state);
    } else {
      // Set state variable without the updated pose
      joint_state_pos_done = joint_state;
      link_state_pos_done = link_state;
    }

    // Init link and joint state for twist and accel computation
    JointState joint_state_twist_done;
    LinkState link_state_twist_done;
    if (compute_twist) {
      joint_state_twist_done =
          forwardJointTwist(joint, parent_link_state_fk_done, joint_state_pos_done);
      link_state_twist_done = forwardLinkTwist(link, joint_state_twist_done, link_state_pos_done);
    } else {
      // Set state variable without the updated twist
      joint_state_twist_done = joint_state_pos_done;
      link_state_twist_done = link_state_pos_done;
    }

    JointState joint_state_accel_done;
    LinkState link_state_accel_done;
    if (compute_accel) {
      joint_state_accel_done =
          forwardJointAccel(joint, parent_link_state_fk_done, joint_state_twist_done);
      link_state_accel_done = forwardLinkAccel(link, joint_state_accel_done, link_state_twist_done);
    } else {
      // Set state variable without the updated accel
      joint_state_accel_done = joint_state_twist_done;
      link_state_accel_done = link_state_twist_done;
    }
    state_variable_new.setJointState(joint_id, joint_state_accel_done);
    state_variable_new.setLinkState(link_id, link_state_accel_done);
  }
  return state_variable_new;
}

Eigen::MatrixXd Kinematics::computeGeneralizedJacobianForLink(const Robot &robot,
                                                              const int link_id) {
  // TODO: This function is not tested enough
  const int DOF = 6;

  auto HH = robot.computeInertiaMatrix();              // Inertia matrix (6+n)x(6+n)
  auto Jb = robot.computeBaseToLinkJacobian(link_id);  // Base to link jacobian
  auto Jm = robot.computeJointToLinkJacobian(link_id); // Joint to link jacobian

  Eigen::MatrixXd Hb = HH.topLeftCorner(DOF, DOF); // Base inertia matrix
  Eigen::MatrixXd Hbm =
      HH.topRightCorner(DOF, robot.getActuatorNumber()); // Base-joint inertia matrix

  Eigen::MatrixXd GJ = Jm - Jb * Hb.inverse() * Hbm; // Generalized Jacobian
  return GJ;
}

Eigen::MatrixXd Kinematics::computeGeneralizedJacobianForEndEffector(const Robot &robot,
                                                                     const int end_effector_id) {
  // TODO: This function is not tested enough
  const int link_id = robot.getEndEffector(end_effector_id).getId();
  return computeGeneralizedJacobianForLink(robot, link_id);
}

Eigen::MatrixXd Kinematics::computeGeneralizedJacobianForEndTip(const Robot &robot,
                                                                const int end_effector_id) {
  // TODO: This function is not tested enough
  const int DOF = 6;

  auto HH = robot.computeInertiaMatrix();                        // Inertia matrix (6+n)x(6+n)
  auto Jb = robot.computeBaseToEndTipJacobian(end_effector_id);  // Base to link jacobian
  auto Jm = robot.computeJointToEndTipJacobian(end_effector_id); // Joint to link jacobian

  Eigen::MatrixXd Hb = HH.topLeftCorner(DOF, DOF); // Base inertia matrix
  Eigen::MatrixXd Hbm =
      HH.topRightCorner(DOF, robot.getActuatorNumber()); // Base-joint inertia matrix

  Eigen::MatrixXd GJ = Jm - Jb * Hb.inverse() * Hbm; // Generalized Jacobian
  return GJ;
}

Eigen::MatrixXd Kinematics::computeJointToLinkJacobian(const Robot &robot, const int link_id) {
  const int DOF = 6;
  const int actuator_number = robot.getActuatorNumber();

  Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(DOF, actuator_number);

  Pose link_pose = robot.getLinkPoseInWorldFrame(link_id);

  Link joint_finder = robot.getLink(link_id); // Used to find parent joint. After the implementation
                                              // of Connection class, this will be removed >> TODO

  int loop_cont = 0; // For safety
  while (joint_finder.getId() > Link::kBase) {
    int joint_id = joint_finder.getParentJointId();
    Joint joint = robot.getJoint(joint_id);
    if (!joint.isActuator()) {
      // Skip fixed joint
      joint_finder = robot.getLink(joint_finder.getParentId());
      continue;
    }
    int actuator_id = joint.getActuatorId();
    JointState joint_state = robot.getJointState(joint_id);

    Eigen::Vector6d jacob_col = Eigen::Vector6d::Zero();

    Eigen::Vector3d trans_joint_to_link =
        joint_state.getPoseInWorldFrame().computeTranslationToPoint(link_pose);

    // Overwrite joint velocity to see its effect on link twist
    // This is the definition of jacobian
    joint_state.setVelocity(1);
    jacob_col =
        joint.twistByActuation(joint_state).computePointTwist(trans_joint_to_link).getTwist();

    jacobian.block(0, actuator_id, DOF, 1) = jacob_col;

    // Move link to its parent
    joint_finder = robot.getLink(joint_finder.getParentId());

    //  Safety to avoid infinite loop
    loop_cont++;
    if (loop_cont > actuator_number) {
      throw std::logic_error("Error: Unexpected error. Loop count is over joint number");
    }
  }
  return jacobian;
}

Eigen::MatrixXd Kinematics::computeJointToEndTipJacobian(const Robot &robot,
                                                         const int end_effector_id) {
  const int DOF = 6;
  const int actuator_number = robot.getActuatorNumber();

  Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(DOF, actuator_number);

  Pose tip_pose = robot.getEndTipPoseInWorldFrame(end_effector_id);

  Link joint_finder =
      robot.getEndEffector(end_effector_id); // Used to find parent joint. After the implementation
                                             // of Connection class, this will be removed >> TODO

  int loop_cont = 0; // For safety
  while (joint_finder.getId() > Link::kBase) {
    int joint_id = joint_finder.getParentJointId();
    Joint joint = robot.getJoint(joint_id);
    if (!joint.isActuator()) {
      // Skip fixed joint
      joint_finder = robot.getLink(joint_finder.getParentId());
      continue;
    }
    int actuator_id = joint.getActuatorId();
    JointState joint_state = robot.getJointState(joint_id);

    Eigen::Vector6d jacob_col = Eigen::Vector6d::Zero();

    Eigen::Vector3d trans_joint_to_link =
        joint_state.getPoseInWorldFrame().computeTranslationToPoint(tip_pose);

    // Overwrite joint velocity to see its effect on end-tip twist
    // This is the definition of jacobian
    joint_state.setVelocity(1);
    jacob_col =
        joint.twistByActuation(joint_state).computePointTwist(trans_joint_to_link).getTwist();

    jacobian.block(0, actuator_id, DOF, 1) = jacob_col;

    // Move link to its parent
    joint_finder = robot.getLink(joint_finder.getParentId());

    //  Safety to avoid infinite loop
    loop_cont++;
    if (loop_cont > actuator_number) {
      throw std::logic_error("Error: Unexpected error. Loop count is over joint number");
    }
  }
  return jacobian;
}

Eigen::MatrixXd Kinematics::computeJointToLinkJacobianDerivative(const Robot &robot,
                                                                 const int link_id) {
  const int DOF = 6;
  const int actuators_number = robot.getActuatorNumber();
  LinkState link_state = robot.getLinkState(link_id);

  Eigen::MatrixXd jacobian_derivative = Eigen::MatrixXd::Zero(DOF, actuators_number);

  Link joint_finder = robot.getLink(link_id); // Used to find parent joint. After the implementation
                                              // of Connection class, this will be removed >> TODO

  int loop_cont = 0; // For safety
  while (joint_finder.getId() > Link::kBase) {
    int joint_id = joint_finder.getParentJointId();
    Joint joint = robot.getJoint(joint_id);
    if (!joint.isActuator()) {
      // Skip fixed joint
      joint_finder = robot.getLink(joint_finder.getParentId());
      continue;
    }
    int actuator_id = joint.getActuatorId();
    JointState joint_state = robot.getJointState(joint_id);

    jacobian_derivative.block(0, actuator_id, DOF, 1) =
        joint.computeVelocityContributionToLinkAccel(joint_state, link_state);

    // Move link to its parent
    joint_finder = robot.getLink(joint_finder.getParentId());

    //  Safety to avoid infinite loop
    loop_cont++;
    if (loop_cont > actuators_number) {
      throw std::logic_error("Error: Unexpected error. Loop count is over joint number");
    }
  }
  return jacobian_derivative;
}

Eigen::MatrixXd Kinematics::computeBaseToLinkJacobian(const Robot &robot, const int link_id) {
  const int DOF = 6;
  Pose link_pose = robot.getLinkPoseInWorldFrame(link_id);
  Pose base_pose = robot.getLinkPoseInWorldFrame(Link::kBase);

  Eigen::MatrixXd jacobian = Eigen::MatrixXd::Identity(DOF, DOF);

  for (int i = 0; i < DOF; i++) {
    Eigen::Vector3d trans_base_to_link = base_pose.computeTranslationToPoint(link_pose);
    jacobian.topRightCorner(3, 3) = -skewSymmetric(trans_base_to_link);
  }

  return jacobian;
}

Eigen::MatrixXd Kinematics::computeBaseToEndTipJacobian(const Robot &robot,
                                                        const int end_effector_id) {
  const int DOF = 6;
  Pose end_tip_pose = robot.getEndTipPoseInWorldFrame(end_effector_id);
  Pose base_pose = robot.getLinkPoseInWorldFrame(Link::kBase);

  Eigen::MatrixXd jacobian = Eigen::MatrixXd::Identity(DOF, DOF);

  for (int i = 0; i < DOF; i++) {
    Eigen::Vector3d trans_base_to_link = base_pose.computeTranslationToPoint(end_tip_pose);
    jacobian.topRightCorner(3, 3) = -skewSymmetric(trans_base_to_link);
  }

  return jacobian;
}

Eigen::MatrixXd Kinematics::computeBaseToLinkJacobianDerivative(const Robot &robot,
                                                                const int link_id) {
  const int DOF = 6;
  Twist link_twist = robot.getLinkState(link_id).getTwistInWorldFrame();
  Twist base_twist = robot.getLinkState(Link::kBase).getTwistInWorldFrame();

  Eigen::MatrixXd jacobian_derivative = Eigen::MatrixXd::Zero(DOF, DOF);

  Eigen::Vector3d relative_velocity =
      link_twist.getLinearVelocity() - base_twist.getLinearVelocity();
  jacobian_derivative.topRightCorner(3, 3) = -skewSymmetric(relative_velocity);

  return jacobian_derivative;
}
} // namespace spacedyn_ros
