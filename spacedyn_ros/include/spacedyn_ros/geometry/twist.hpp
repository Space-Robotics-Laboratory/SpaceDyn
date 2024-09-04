#ifndef SPACEDYN_ROS_TWIST_HPP_
#define SPACEDYN_ROS_TWIST_HPP_

#include "eigen3/Eigen/Core"
#include "geometry_msgs/msg/twist.hpp"
#include "spacedyn_ros/geometry/pose.hpp"

namespace spacedyn_ros {
class Twist {
private:
  Eigen::VectorXd twist_;
  Frame frame_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Twist(const Frame &frame = Frame::kWorld,
        const Eigen::VectorXd &twist = Eigen::VectorXd::Zero(6));
  Twist(const Frame &frame, const Eigen::Vector3d &linier_velocity,
        const Eigen::Vector3d &angular_velocity);
  ~Twist() = default;

  const Frame &getFrame() const;

  /**
   * @fn getOriginTwist()
   * @brief Get the twist of the link in World frame. Do not confuse it with the
   * twist used in the Screw theory.
   */
  const Eigen::VectorXd &getOriginTwist() const;
  /**
   * @fn getVelocity()
   * @brief Get the linear velocity of the link in World frame
   */
  Eigen::Vector3d getOriginLinierVelocity() const;
  /**
   * @fn getAngularVelocity()
   * @brief Get the angular velocity vector of the link in World frame.
   */
  Eigen::Vector3d getOriginAngularVelocity() const;

  Twist operator+(const Twist &twist) const;
  Twist operator-(const Twist &twist) const;

  /**
   * @fn computePointTwist()
   * @brief Get the twist of the point by this object in World frame:
   * [vs, ws] =
   * [v + w x r, w]
   * @param tf_to_point: the position of the point in the world frame
   * @return the twist of the point in the world frame
   */
  Twist computePointTwist(const Transform &tf_to_point) const; // This function is not complete
  // TODO: Remove this function and use Twist::computePointTwist(const Transform &tf_to_point)
  Twist computePointTwist(const Eigen::Vector3d &translation_to_point_in_world_frame) const;

  // ROS Interface
  geometry_msgs::msg::Twist toRosMessage() const;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_TWIST_HPP_
