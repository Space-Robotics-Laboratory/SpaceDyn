#ifndef SPACEDYN_ROS_POSE_HPP_
#define SPACEDYN_ROS_POSE_HPP_

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"
#include "geometry_msgs/msg/pose.hpp"
#include "geometry_msgs/msg/transform_stamped.hpp"
#include "spacedyn_ros/geometry/frame.hpp"
#include "spacedyn_ros/geometry/transform.hpp"
#include "spacedyn_ros/geometry/twist.hpp"
#include "spacedyn_ros/geometry/wrench.hpp"

namespace spacedyn_ros {
class Twist;
class Pose {
private:
  Eigen::Isometry3d pose_in_world_;
  void checkRotationMatrixNormalized(const Eigen::Matrix3d &rotation_matrix) const;
  void checkQuaternionNormalized(const Eigen::Quaterniond &quaternion) const;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Pose(const Eigen::Isometry3d &pose = Eigen::Isometry3d::Identity());
  Pose(const Eigen::Vector3d &position, const Eigen::Matrix3d &attitude);
  Pose(const Eigen::Vector3d &position, const Eigen::Quaterniond &attitude);
  // ROS Interface
  Pose(const geometry_msgs::msg::Pose &pose_msg);
  ~Pose() = default;

  /**
   * @fn getPose()
   * @brief Get the pose of the link in World frame
   */
  const Eigen::Isometry3d &getOriginPose() const;

  /**
   * @fn getOriginPosition()
   * @brief Get the position of the link in World frame
   */
  Eigen::Vector3d getOriginPosition() const;

  /**
   * @fn getOriginAttitude()
   * @brief Get the attitude of the link in World frame
   */
  Eigen::MatrixXd getOriginAttitude() const;
  Eigen::Quaterniond getOriginQuaternion() const;

  Eigen::Vector3d computeTranslationToPoint(const Pose &point_pose) const;
  Transform computeTransformToPoint(const Frame &frame, const Pose &point_pose) const;
  Pose computePointPose(const Transform &tf_to_point) const;
  Pose computeOriginPoseFromPointPose(const Transform &tf_to_point, const Pose &point_pose) const;

  Twist computeTwistInLocalFrame(const Twist &twist_in_world_frame) const;
  Wrench computeWrenchInWorldFrame(const Wrench &wrench_in_local_frame) const;

  Eigen::Quaterniond computeDerivativeAttitude(const Twist &twist) const;
  // ROS Interface
  geometry_msgs::msg::Pose toRosMessage() const;
  geometry_msgs::msg::TransformStamped toRosMessage(const std::string &frame_name,
                                                    const std::string &child_frame_name) const;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_POSE_HPP_
