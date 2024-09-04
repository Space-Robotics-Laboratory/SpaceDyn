#ifndef SPACEDYN_ROS_ACCEL_HPP_
#define SPACEDYN_ROS_ACCEL_HPP_

#include "eigen3/Eigen/Core"
#include "geometry_msgs/msg/accel.hpp"
#include "spacedyn_ros/geometry/twist.hpp"

namespace spacedyn_ros {
class Accel {
private:
  Eigen::VectorXd accel_;
  Frame frame_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Accel(const Frame &frame = Frame::kWorld,
        const Eigen::VectorXd &accel = Eigen::VectorXd::Zero(6));
  ~Accel() = default;

  const Frame &getFrame() const;
  const Eigen::VectorXd &getOriginAccel() const;
  Eigen::Vector3d getOriginLinierAcceleration() const;
  Eigen::Vector3d getOriginAngularAcceleration() const;

  Accel operator+(const Accel &accel) const;
  Accel operator-(const Accel &accel) const;

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
