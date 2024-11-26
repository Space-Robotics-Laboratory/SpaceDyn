#include "spd_example/zhanibekov_effect_sim.hpp"
#include "ament_index_cpp/get_package_share_directory.hpp"
#include "spacedyn_ros/linkage/linkage.hpp"
#include <rclcpp/rclcpp.hpp>

namespace spacedyn_ros {
ZhanibekovEffectSimulation::ZhanibekovEffectSimulation(const std::string &name,
                                                       const std::string &path_to_urdf)
    : Node(name) {
  timer_ = this->create_wall_timer(std::chrono::milliseconds(1),
                                   std::bind(&ZhanibekovEffectSimulation::timer_callback, this));
  tf_broadcaster_ = std::make_unique<tf2_ros::TransformBroadcaster>(this);

  Model model(path_to_urdf);
  model.setGravity(Eigen::Vector3d(0, 0, 0));
  model.setDeltaTimeMicroSec(1000);
  this->robot_ = Robot(model);

  // Apply external wrench to the base link
  Eigen::VectorXd wrench = Eigen::VectorXd::Zero(6);
  wrench.tail(3) = Eigen::Vector3d(1, 5000, 0);
  robot_.applyExternalWrench(Link::ID::kBase, Wrench(Frame::kWorld, wrench));
}

void ZhanibekovEffectSimulation::timer_callback() {
  robot_.step();
  auto tf_base = robot_.basePoseToRosTf();
  tf_broadcaster_->sendTransform(tf_base);
  auto tf_joints = robot_.jointPoseToRosTf();
  tf_base.header.stamp = this->now();
  for (auto &tf_joint : tf_joints) {
    tf_joint.header.stamp = this->now();
    tf_broadcaster_->sendTransform(tf_joint);
  }
}
} // namespace spacedyn_ros

int main(int argc, char *argv[]) {
  rclcpp::init(argc, argv);
  std::string path_to_urdf = ament_index_cpp::get_package_share_directory("spd_example") +
                             "/urdf/zhanibekov_effect_model.urdf";
  auto node = std::make_shared<spacedyn_ros::ZhanibekovEffectSimulation>("zhanibekov_effect_sim",
                                                                         path_to_urdf);
  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}
