#ifndef SPACE_DYN_CONNECTION_HPP_
#define SPACE_DYN_CONNECTION_HPP_
#include "eigen3/Eigen/Core"

namespace spacedyn_ros {
class Connection {
private:
  Eigen::Sparse connection_jag_;

public:
  Connection(/* args */) {}
  ~Connection() = default;
};
} // namespace spacedyn_ros

#endif // SPACE_DYN_CONNECTION_HPP_