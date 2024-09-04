#ifndef SPACEDYN_ROS_INERTIA_HPP_
#define SPACEDYN_ROS_INERTIA_HPP_

#include "eigen3/Eigen/Core"
#include "spacedyn_ros/geometry/accel.hpp"
#include "spacedyn_ros/geometry/frame.hpp"
#include "spacedyn_ros/geometry/wrench.hpp"

#include "geometry_msgs/msg/inertia.hpp"

namespace spacedyn_ros {
class Inertia {
private:
  double mass_;
  Eigen::Matrix3d inertia_;
  Frame frame_;
  void checkMass(const double mass) const;
  void checkInertia(const Eigen::Matrix3d &inertia) const;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Inertia(const Frame &frame = Frame::kLocal, const double mass = 1,
          const Eigen::Matrix3d &inertia = Eigen::Matrix3d::Identity(3, 3));
  ~Inertia() = default;

  double getMass() const;
  const Eigen::Matrix3d &getOriginInertiaTensor() const;
  const Frame &getFrame() const;

  // ROS Interface
  geometry_msgs::msg::Inertia toRosMessage() const;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_INERTIA_HPP_
