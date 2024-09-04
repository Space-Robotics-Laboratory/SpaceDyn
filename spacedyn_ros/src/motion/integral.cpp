#include "spacedyn_ros/motion/integral.hpp"
#include "spacedyn_ros/motion/dynamics.hpp"
#include "spacedyn_ros/motion/kinematics.hpp"
#include "spacedyn_ros/robot/robot.hpp"
#include <iostream>

namespace spacedyn_ros {
Integral::Integral() {}

Eigen::Quaterniond Integral::addScaledQuaternion(const double scale_a, const Eigen::Quaterniond &a,
                                                 const double scale_b,
                                                 const Eigen::Quaterniond &b) {
  Eigen::Quaterniond result =
      Eigen::Quaterniond(scale_a * a.w() + scale_b * b.w(), scale_a * a.x() + scale_b * b.x(),
                         scale_a * a.y() + scale_b * b.y(), scale_a * a.z() + scale_b * b.z())
          .normalized();
  return result;
}

StateVariable Integral::stepGeneralizedState(const StateVariable &sv, const StateVariable &dsv,
                                             const double dt) {
  StateVariable state_variable_next = sv;
  Pose base_pose = sv.getLinkState(Link::ID::kBase).getPoseInWorldFrame();
  Twist base_twist = dsv.getLinkState(Link::ID::kBase).getTwistInWorldFrame();

  // Compute the next base pose and joint position
  auto base_quat = base_pose.getOriginQuaternion();
  auto base_quat_d = base_pose.computeDerivativeAttitude(base_twist);
  Eigen::Vector3d base_pos_next =
      base_pose.getOriginPosition() + base_twist.getOriginLinierVelocity() * dt;
  Eigen::Quaterniond base_quat_next = addScaledQuaternion(1.0, base_quat, dt, base_quat_d);
  Eigen::Isometry3d base_pose_next = Eigen::Translation3d(base_pos_next) * base_quat_next;
  Eigen::VectorXd joint_pose_next = sv.getJointPosition() + dsv.getJointVelocity() * dt;
  state_variable_next.setGeneralizedCoordinates(base_pose_next, joint_pose_next);
  state_variable_next.setGeneralizedVelocity(sv.getGeneralizedVelocity() +
                                             dsv.getGeneralizedAcceleration() * dt);
  return state_variable_next;
}

StateVariable Integral::addScaledGeneralizedVelocity(const double scale_a, const StateVariable &a,
                                                     const double scale_b, const StateVariable &b) {
  StateVariable result = a;
  Eigen::VectorXd vel_a = a.getGeneralizedVelocity();
  Eigen::VectorXd vel_b = b.getGeneralizedVelocity();
  result.setGeneralizedVelocity(vel_a * scale_a + vel_b * scale_b);
  return result;
}

StateVariable Integral::addScaledGeneralizedAcceleration(const double scale_a,
                                                         const StateVariable &a,
                                                         const double scale_b,
                                                         const StateVariable &b) {
  StateVariable result = a;
  Eigen::VectorXd accel_a = a.getGeneralizedAcceleration();
  Eigen::VectorXd accel_b = b.getGeneralizedAcceleration();
  result.setGeneralizedAcceleration(accel_a * scale_a + accel_b * scale_b);
  return result;
}

StateVariable Integral::sumUpWeighedDerivativeOfGeneralizedStates(
    const std::vector<std::pair<double, StateVariable>> &args) {
  StateVariable result = args[0].second.copyVacant();
  for (size_t i = 0; i < args.size(); i++) {
    result = addScaledGeneralizedVelocity(1, result, args[i].first, args[i].second);     // vel
    result = addScaledGeneralizedAcceleration(1, result, args[i].first, args[i].second); // accel
  }
  return result;
}

StateVariable Integral::rungeKutta4(const Robot &robot) {
  const double dt = robot.getModel().getDeltaTimeSec();
  const double half_dt = dt / 2.0;
  const double sixth_dt = dt / 6.0;

  // Only derivative of the generalized states are considered
  const StateVariable sv0 = robot.computeForward(); // All states are up2date to t = t0

  // k1 = f(sv0)
  const StateVariable k1 = sv0;
  // k2 = f(sv0 + k1 * dt / 2)
  const StateVariable k2 = Robot(robot, stepGeneralizedState(sv0, k1, half_dt)).computeForward();
  // k3 = f(sv0 + k2 * dt / 2)
  const StateVariable k3 = Robot(robot, stepGeneralizedState(sv0, k2, half_dt)).computeForward();
  // k4 = f(sv0 + k3 * dt)
  const StateVariable k4 = Robot(robot, stepGeneralizedState(sv0, k3, dt)).computeForward();

  // sv_next = sv0 + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
  const std::vector<std::pair<double, StateVariable>> args = {{1, k1}, {2, k2}, {2, k3}, {1, k4}};
  const StateVariable dsv = sumUpWeighedDerivativeOfGeneralizedStates(args);

  const StateVariable sv_new = stepGeneralizedState(sv0, dsv, sixth_dt);

  return Kinematics::computeForward(Robot(robot, sv_new), true, true, false);
}

StateVariable Integral::euler(const Robot &robot) {
  const double dt = robot.getModel().getDeltaTimeSec();
  const StateVariable sv_0 = robot.computeForward(); // First update accel
  const StateVariable sv_new = stepGeneralizedState(sv_0, sv_0, dt);
  return Kinematics::computeForward(Robot(robot, sv_new), true, true, false);
}

} // namespace spacedyn_ros
