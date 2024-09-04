#ifndef SPACEDYN_ROS_MODEL_HPP_
#define SPACEDYN_ROS_MODEL_HPP_

#include "eigen3/Eigen/Core"
#include "spacedyn_ros/linkage/linkage.hpp"

namespace spacedyn_ros {
class Model {
private:
  size_t dt_microsec_;
  Eigen::Vector3d gravity_;
  Linkage linkage_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Model(const Linkage &linkage = Linkage());
  Model(const std::string &path_to_urdf);
  ~Model() = default;

  void loadURDF(const std::string &path_to_urdf);

  /**
   * @fn
   * Write the function description here
   * @brief Summary description
   * @param (parameter name) Description of the parameter
   * @param (parameter name) Description of the parameter
   * @return Description of the return value
   * @sa If you write the function to refer to, a link can be created
   * @detail Detailed description
   */
  const Linkage &getLinkage() const;
  double getTotalMass() const;
  int getLinkNumber() const;
  int getJointNumber() const;

  double getDeltaTimeSec() const;
  double getDeltaTimeMilliSec() const;
  size_t getDeltaTimeMicroSec() const;

  const Eigen::Vector3d &getGravity() const;

  void setDeltaTimeMicroSec(const size_t dt_microsec);
  void setGravity(const Eigen::Vector3d &gravity);
};
} // namespace spacedyn_ros

#endif // SPACEDYN_ROS_MODEL_HPP_
