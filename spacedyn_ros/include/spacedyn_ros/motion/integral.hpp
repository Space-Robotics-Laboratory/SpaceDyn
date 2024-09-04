#ifndef SPACEDYN_ROS_INTEGRAL_HPP_
#define SPACEDYN_ROS_INTEGRAL_HPP_
#include "spacedyn_ros/robot/robot.hpp"
#include "spacedyn_ros/robot/state_variable.hpp"

namespace spacedyn_ros {
class Integral {
private:
  /**
   * @brief Add two quaternions. Used in stepGeneralizedState method.
   */
  static Eigen::Quaterniond addScaledQuaternion(const double scale_a, const Eigen::Quaterniond &a,
                                                const double scale_b, const Eigen::Quaterniond &b);

  /**
   * @brief Compute the next state variable from the current state variable. Only the generalized
   * coordinates and velocities are updated. sv_next = sv + dsv * dt
   */
  static StateVariable stepGeneralizedState(const StateVariable &sv, const StateVariable &dsv,
                                            const double dt);

  /**
   * @brief Add two StateVariable objects. Only the addGeneralizedVelocity states are added. Used in
   * Runge-Kutta 4th order method.
   */
  static StateVariable addScaledGeneralizedVelocity(const double scale_a, const StateVariable &a,
                                                    const double scale_b, const StateVariable &b);

  /**
   * @brief Add two StateVariable objects. Only the addGeneralizedAcceleration states are added.
   * Used in Runge-Kutta 4th order method.
   */
  static StateVariable addScaledGeneralizedAcceleration(const double scale_a,
                                                        const StateVariable &a,
                                                        const double scale_b,
                                                        const StateVariable &b);

  /**
   * @brief Sum up the weighted derivative of generalized states. args = { {weight, dsv1}, ... }.
   * Used in Runge-Kutta 4th order method.
   */
  static StateVariable sumUpWeighedDerivativeOfGeneralizedStates(
      const std::vector<std::pair<double, StateVariable>> &args);

public:
  Integral();
  ~Integral() = default;
  static StateVariable rungeKutta4(const Robot &robot);
  static StateVariable euler(const Robot &robot);
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_INTEGRAL_HPP_
