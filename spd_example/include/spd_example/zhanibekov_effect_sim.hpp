#ifndef SPACEDYN_ROS_CUBIC_SIMULATION_HPP_
#define SPACEDYN_ROS_CUBIC_SIMULATION_HPP_
#include "rclcpp/rclcpp.hpp"
#include "spacedyn_ros/robot/robot.hpp"
#include "tf2_ros/transform_broadcaster.h"

namespace spacedyn_ros {
class ZhanibekovEffectSimulation : public rclcpp::Node {
public:
  ZhanibekovEffectSimulation(const std::string &name, const std::string &path_to_urdf);
  ~ZhanibekovEffectSimulation() = default;

private:
  Robot robot_;
  rclcpp::TimerBase::SharedPtr timer_;
  void timer_callback();

  std::unique_ptr<tf2_ros::TransformBroadcaster> tf_broadcaster_;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_CUBIC_SIMULATION_HPP_
