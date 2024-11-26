#ifndef SPACE_DYN_ROS_CONNECTION_HPP_
#define SPACE_DYN_ROS_CONNECTION_HPP_

#include "spacedyn_ros/linkage/link.hpp"
#include <boost/assign/list_of.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <deque>
#include <eigen3/Eigen/Core>
#include <iostream>
#include <map>

namespace spacedyn_ros {
class Connection {
public:
  Connection();
  ~Connection() = default;

  // Define graph properties
  typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, boost::no_property,
                                boost::property<boost::edge_weight_t, int>>
      Graph;
  Graph link_connection_graph; // express connection of robot links

  /**
   * @brief Add node as base to link_connection_graph
   *
   */
  void addBase();

  /**
   * @brief Add node as link to link_connection_graph
   *
   * @param parent_link_id
   * @param child_link_id
   */
  void addLink(const int parent_link_id, const int child_link_id);

  /**
   * @brief Get the Link Id Chain object
   *
   * @param graph
   * @param start_link_id
   * @param end_link_id
   * @return std::vector<int> id_chain
   */
  std::vector<int> getLinkIdChain(const int start_link_id, const int end_link_id) const;

private:
  /**
   * @brief Check input link IDs are valid or not
   *
   * @param parent_link_id
   * @param child_link_id
   * @return true
   * @return false
   */
  bool isValidNewLinkId(const int parent_link_id, const int child_link_id) const;

  /**
   * @brief Check inputted link ID exists in link connection graph
   *
   * @param link_id
   * @return true
   * @return false
   */
  bool linkIdExists(const int link_id) const;

  /**
   * @brief Get Reversed Map object for link ID and vertex ID
   *
   * @param map
   * @return std::map<Graph::vertex_descriptor, int>
   */
  std::map<Graph::vertex_descriptor, int>
  getReverseMap(const std::map<int, Graph::vertex_descriptor> &map) const;

  /**
   * @brief Visualize properties of graph
   *
   * @param Graph g
   */
  void visualizeGraph(const Graph &g) const;

  std::map<int, Graph::vertex_descriptor> v_desc; // contain link ID and vertex ID
};
} // namespace spacedyn_ros

#endif // SPACE_DYN_ROS_CONNECTION_HPP_
