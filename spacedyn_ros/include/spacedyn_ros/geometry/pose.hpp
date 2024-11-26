#ifndef SPACEDYN_ROS_POSE_HPP_
#define SPACEDYN_ROS_POSE_HPP_

#include "spacedyn_ros/geometry/frame.hpp"
#include "spacedyn_ros/geometry/transform.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <geometry_msgs/msg/pose.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>

namespace spacedyn_ros {
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
  ~Pose() = default;

  /**
   * @fn getPoseInWorldFrame()
   * @brief Get the pose of the link in World frame
   */
  const Eigen::Isometry3d &getPoseInWorldFrame() const;

  /**
   * @fn getPositionInWorldFrame()
   * @brief Get the position of the link in World frame
   */
  Eigen::Vector3d getPositionInWorldFrame() const;

  /**
   * @fn getAttitudeInWorldFrame()
   * @brief Get the attitude of the link in World frame
   */
  Eigen::Matrix3d getAttitudeInWorldFrame() const;
  Eigen::Quaterniond getQuaternionInWorldFrame() const;

  Eigen::Vector3d computeTranslationToPoint(const Pose &point_pose) const;
  Transform computeTransformToPoint(const Frame &frame, const Pose &point_pose) const;
  Pose computePointPose(const Transform &tf_to_point) const;
  Pose computePoseByInvertingPointPose(const Transform &tf_to_point, const Pose &point_pose) const;

  // ROS Interface
  geometry_msgs::msg::Pose toRosMessage() const;
  geometry_msgs::msg::TransformStamped toRosMessage(const std::string &frame_name,
                                                    const std::string &child_frame_name) const;
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_POSE_HPP_
