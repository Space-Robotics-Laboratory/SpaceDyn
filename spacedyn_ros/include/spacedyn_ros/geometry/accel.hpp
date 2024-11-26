#ifndef SPACEDYN_ROS_ACCEL_HPP_
#define SPACEDYN_ROS_ACCEL_HPP_

#include "spacedyn_ros/geometry/twist.hpp"
#include "spacedyn_ros/util/matrix_operation.hpp"
#include <eigen3/Eigen/Core>
#include <geometry_msgs/msg/accel.hpp>

namespace spacedyn_ros {
class Accel {
private:
  Eigen::Vector6d accel_;
  Frame frame_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Accel(const Frame &frame = Frame::kWorld,
        const Eigen::Vector6d &accel = Eigen::Vector6d::Zero(6));
  Accel(const Frame &frame, const Eigen::Vector3d &linear_acceleration,
        const Eigen::Vector3d &angular_acceleration);
  ~Accel() = default;

  const Frame &getFrame() const;
  const Eigen::Vector6d &getAccel() const;
  Eigen::Vector3d getLinearAcceleration() const;
  Eigen::Vector3d getAngularAcceleration() const;

  Accel operator+(const Accel &accel) const;
  Accel operator-(const Accel &accel) const;
  Accel getAccelInFrame(const Frame &frame, const Pose &pose) const;

  // TODO: Implement the following functions
  Accel computePointAccel(const Twist &twist, const Transform &tf_to_point) const;

  // TODO: Stop using Vector3d and use Transform instead
  Accel computePointAccel(const Twist &twist,
                          const Eigen::Vector3d &translation_to_point_in_world_frame) const;

  // ROS Interface
  geometry_msgs::msg::Accel toRosMessage() const;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_ACCEL_HPP_
