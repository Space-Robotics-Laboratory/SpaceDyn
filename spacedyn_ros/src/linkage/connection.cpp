#include "spacedyn_ros/linkage/connection.hpp"

#define DEBUG false

namespace spacedyn_ros {
Connection::Connection(/* args */) {}

void Connection::visualizeGraph(const Graph &g) const {
  std::cout << "visualize graph" << std::endl;

  // Output vertex information
  auto v = boost::vertices(g);
  for (auto i = v.first; i != v.second; i++) {
    std::cout << "vertex: " << *i << std::endl;
  }
  std::cout << "num_vertices = " << boost::num_vertices(g) << std::endl;

  // Visualize graph
#if DEBUG
  const std::string file_name = "boost_graph.dot";
  std::ofstream file(file_name);
  if (!file) {
    std::cout << "Opening file failed." << std::endl;
    return;
  }
  boost::write_graphviz(file, g);
#endif
}

bool Connection::isValidNewLinkId(const int parent_link_id, const int child_link_id) const {
  if (parent_link_id == child_link_id) {
    throw std::invalid_argument("Error: Parent link ID is " + std::to_string(parent_link_id) +
                                " and child link ID is also " + std::to_string(child_link_id));
  } else if (boost::num_vertices(link_connection_graph) == 0) {
    throw std::runtime_error("Error: Base link is necessary to add child link.");
  } else if (!linkIdExists(parent_link_id)) {
    throw std::runtime_error("Error: parent_link_id does not exist.");
  } else {
    // TODO: Below part should be deleted after closed loop is implemented.
    bool is_unique_child_link_id = true;
    for (auto i = v_desc.begin(); i != v_desc.end(); i++) {
      if (child_link_id == i->first) {
        is_unique_child_link_id = false;
      }
    }
    if (!is_unique_child_link_id) {
      throw std::runtime_error("Error: Closed loop link is not allowed.");
    }
  }
  return true;
}

void Connection::addBase() {
  // Base link must have link ID as 0.
  if (boost::num_vertices(link_connection_graph) != 0) {
    throw std::runtime_error("Error: Base ID must be 0.");
  }

  const int base_id = 0;
  v_desc.insert(std::make_pair(base_id, add_vertex(link_connection_graph)));
}

void Connection::addLink(const int parent_link_id, const int child_link_id) {
  try {
    isValidNewLinkId(parent_link_id, child_link_id);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    throw std::runtime_error("Error: Invalid ID is found.");
  }

  // Add vertex of child link.
  v_desc.insert(std::make_pair(child_link_id, add_vertex(link_connection_graph)));

  // Add edge between parent and child link.
  const int weight = 1; // All weight is 1 for now.
  add_edge(v_desc.at(parent_link_id), v_desc.at(child_link_id), weight, link_connection_graph);

#if DEBUG
  std::cout << "parent_link_id = " << parent_link_id << std::endl;
  std::cout << "child_link_id = " << child_link_id << std::endl;
  std::cout << "v_desc.at(" << parent_link_id << ") = " << v_desc.at(parent_link_id) << std::endl;
  std::cout << "v_desc.at(" << child_link_id << ") = " << v_desc.at(child_link_id) << std::endl;

  boost::print_graph(link_connection_graph);
#endif
}

std::map<Connection::Graph::vertex_descriptor, int>
Connection::getReverseMap(const std::map<int, Graph::vertex_descriptor> &map) const {
  std::map<Graph::vertex_descriptor, int> rev_map;
  for (auto &pair : map) {
    rev_map[pair.second] = pair.first;
  }
  return rev_map;
}

bool Connection::linkIdExists(const int link_id) const {
  for (auto i = v_desc.begin(); i != v_desc.end(); i++) {
    if (link_id == i->first) {
      return true;
    }
  }
  return false;
}

std::vector<int> Connection::getLinkIdChain(const int start_link_id, const int end_link_id) const {
  if (start_link_id < 0 || end_link_id < 0) {
    throw std::invalid_argument("Error: Inputs must be natural numbers.");
  }

  if (!linkIdExists(start_link_id)) {
    throw std::runtime_error("Error: start_link_id doesn't exist in graph.");
  } else if (!linkIdExists(end_link_id)) {
    throw std::runtime_error("Error: end_link_id doesn't exist in graph.");
  }

  // Containers for distances and predecessors
  std::vector<int> distances(boost::num_vertices(link_connection_graph));
  std::vector<boost::graph_traits<Graph>::vertex_descriptor> predecessors(
      boost::num_vertices(link_connection_graph));

  // Calculate the shortest path
  boost::dijkstra_shortest_paths(
      link_connection_graph, v_desc.at(start_link_id),
      boost::predecessor_map(&predecessors.at(0)).distance_map(&distances.at(0)));

  // Create path from start to end link IDs.
  std::deque<int> path;
  for (int i = v_desc.at(end_link_id); i != v_desc.at(start_link_id); i = predecessors.at(i)) {
    path.push_front(i);
  }
  path.push_front(v_desc.at(start_link_id));

  // Convert std::deque to std::vector
  auto v_desc_reverse = getReverseMap(v_desc);

  std::vector<int> id_chain_from_start_to_end(path.size());
  for (int i = 0; i < path.size(); i++) {
    // link ID and graph node ID are different.
    id_chain_from_start_to_end.at(i) = v_desc_reverse.at(path.at(i));
  }

  return id_chain_from_start_to_end;
}

} // namespace spacedyn_ros
