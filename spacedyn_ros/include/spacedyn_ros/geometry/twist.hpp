#ifndef SPACEDYN_ROS_TWIST_HPP_
#define SPACEDYN_ROS_TWIST_HPP_

#include "spacedyn_ros/geometry/frame.hpp"
#include "spacedyn_ros/geometry/pose.hpp"
#include "spacedyn_ros/geometry/transform.hpp"
#include "spacedyn_ros/util/matrix_operation.hpp"
#include <eigen3/Eigen/Core>
#include <geometry_msgs/msg/twist.hpp>

namespace spacedyn_ros {
class Twist {
private:
  Eigen::Vector6d twist_;
  Frame frame_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Twist(const Frame &frame = Frame::kWorld,
        const Eigen::Vector6d &twist = Eigen::Vector6d::Zero(6));
  Twist(const Frame &frame, const Eigen::Vector3d &linear_velocity,
        const Eigen::Vector3d &angular_velocity);
  ~Twist() = default;

  const Frame &getFrame() const;

  /**
   * @fn getTwist()
   * @brief Get the twist of the link in its frame. Do not confuse it with the
   * twist used in the Screw theory.
   */
  const Eigen::Vector6d &getTwist() const;
  /**
   * @fn getVelocity()
   * @brief Get the linear velocity of the link in its frame
   */
  Eigen::Vector3d getLinearVelocity() const;
  /**
   * @fn getAngularVelocity()
   * @brief Get the angular velocity vector of the link in its frame.
   */
  Eigen::Vector3d getAngularVelocity() const;

  Twist operator+(const Twist &twist) const;
  Twist operator-(const Twist &twist) const;
  Twist getTwistInFrame(const Frame &frame, const Pose &pose) const;

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

  Eigen::Quaterniond computeDerivativeAttitude(const Pose &pose) const;

  // ROS Interface
  geometry_msgs::msg::Twist toRosMessage() const;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_TWIST_HPP_
