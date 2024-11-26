#ifndef SPACEDYN_ROS_WRENCH_HPP_
#define SPACEDYN_ROS_WRENCH_HPP_

#include "spacedyn_ros/geometry/frame.hpp"
#include "spacedyn_ros/geometry/pose.hpp"
#include "spacedyn_ros/geometry/transform.hpp"
#include "spacedyn_ros/util/matrix_operation.hpp"
#include <eigen3/Eigen/Core>
#include <geometry_msgs/msg/wrench.hpp>

namespace spacedyn_ros {
class Wrench {
private:
  Frame frame_;
  Eigen::Vector6d wrench_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Wrench(const Frame &frame = Frame::kWorld,
         const Eigen::Vector6d &wrench = Eigen::Vector6d::Zero(6));
  Wrench(const Frame &frame, const Eigen::Vector3d &force, const Eigen::Vector3d &torque);
  ~Wrench() = default;

  const Frame &getFrame() const;

  const Eigen::Vector6d &getWrench() const;

  Eigen::Vector3d getForce() const;

  Eigen::Vector3d getTorque() const;

  Wrench operator+(const Wrench &wrench) const;
  Wrench operator-() const;
  Wrench operator-(const Wrench &wrench) const;
  Wrench operator+=(const Wrench &wrench);
  Wrench getWrenchInFrame(const Frame &frame, const Pose &pose) const;

  // This function is not complete
  Wrench computePointWrench(const Transform &tf_to_point) const;
  // TODO: Remove this function and use Wrench::computePointWrench(const Transform &tf_to_point)
  /**
   * @fn computeWrenchByInvertingPointWrench
   * @brief Calculate Ta = Pab x Fb
   *
   * @param translation_to_point_in_world_frame
   * @return Wrench
   */
  Wrench computeWrenchByInvertingPointWrench(
      const Eigen::Vector3d &translation_to_point_in_world_frame) const;

  // ROS Interface
  geometry_msgs::msg::Wrench toRosMessage() const;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_WRENCH_HPP_
