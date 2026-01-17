#include <Rcpp.h>
#include <algorithm>  // for std::reverse
using namespace Rcpp;

// [[Rcpp::export]]
LogicalMatrix checkPTA(NumericVector pta_nodes, List points) {

  int k = pta_nodes.size();
  int n = points.size();

  LogicalMatrix out(n, k);

  for (int i = 0; i < n; i++) {
    NumericVector point = points[i];

    // Convert "point" to std::unordered_set for faster look-up
    std::unordered_set<double> set_i(point.begin(), point.end());

    for (int j = 0; j < k; j++) {
      out(i, j) = set_i.count(pta_nodes[j]) > 0;
    }
  }

  CharacterVector colnames = as<CharacterVector>(wrap(pta_nodes));

  CharacterVector rownames = as<CharacterVector>(wrap(points.names()));

  out.attr("dimnames") = List::create(rownames, colnames);


  return out;
}

// [[Rcpp::export]]
List get_node_paths(List rpart_obj) {
  // Extract the frame from the rpart object
  DataFrame frame = as<DataFrame>(rpart_obj["frame"]);

  // Node IDs are stored as row names
  IntegerVector node_ids = frame.attr("row.names");

  int n_nodes = node_ids.size();
  List paths(n_nodes);
  CharacterVector names(n_nodes);

  for (int i = 0; i < n_nodes; ++i) {
    int current = node_ids[i];
    std::vector<int> path;

    // Optional but sensible for rpart trees
    path.reserve(32);

    // Build path bottom-up
    while (current >= 1) {
      path.push_back(current);

      if (current == 1) {
        break;
      }
      current /= 2;
    }

    // Reverse to get root -> node order
    std::reverse(path.begin(), path.end());

    paths[i] = wrap(path);
    names[i] = std::to_string(node_ids[i]);
  }

  paths.attr("names") = names;
  return paths;
}
